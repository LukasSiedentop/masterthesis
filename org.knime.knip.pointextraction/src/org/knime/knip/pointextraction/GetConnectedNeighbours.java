package org.knime.knip.pointextraction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Set;

import org.knime.core.data.def.IntCell;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.view.Views;

@Plugin(menu = {@Menu(label = "PointExtraction"),
                @Menu(label = "Connected Neighbours")}, description = "Gets the Neighbours of Nodes which are connected by a pixel line. Takes a skeleton BitType Image and a Labeling with the Nodes.", headless = true, type = Command.class)
public class GetConnectedNeighbours implements Command {

        @Parameter(type = ItemIO.INPUT, label = "(skeletonized) Image", description = "Skeletonized Image whose Intersections are to be found.")
        private ImgPlus<BitType> input;

        @Parameter(type = ItemIO.INPUT, label = "Labeling", description = "Labeling containing Labels at each intersection (aka node).")
        private ImgLabeling<Integer, ? extends IntegerType<?>> labeling;

        // geht nicht
        //@Parameter(type = ItemIO.INPUT)
        //private ExecutionContext exec;

        //@Parameter(type = ItemIO.OUTPUT)
        //private BufferedDataTable output;

        @Parameter(type = ItemIO.OUTPUT)
        private ImgPlus<BitType> output;

        //private RandomAccess<LabelingType<Integer>> labelingRndAccess = labeling.randomAccess();

        @Override
        public void run() {
                // to eliminate error message
                output = input.copy();

                RandomAccess<LabelingType<Integer>> labelingRndAccess = labeling.randomAccess();

                LabelRegions<Integer> regions = new LabelRegions<>(labeling);
                Set<Integer> nodes = regions.getExistingLabels();

                // 3^n ROI generieren
                final RectangleShape roi = new RectangleShape(1, false);
                RandomAccess<Neighborhood<BitType>> roiAccess = roi.neighborhoodsRandomAccessible(Views.extendValue(input, new BitType(false)))
                                .randomAccess();

                // Positichonsarrays
                int[] nodePixelPosition = new int[labeling.numDimensions()];
                int[] nextOnLine = new int[labeling.numDimensions()];

                // Centerdaten
                double[][] centers = new double[nodes.size() + 1][labeling.numDimensions()];

                // Datenstruktur um Ergebnis aufzunehmen
                @SuppressWarnings("unchecked")
                ArrayList<IntCell>[] nodeArray = ((ArrayList<IntCell>[]) new ArrayList[nodes.size() + 1]);

                // Iteration über alle Knoten
                for (Integer node : nodes) {

                        nodeArray[node] = new ArrayList<IntCell>(0);

                        // Cursor über Knotenpixel
                        Cursor<Void> nodeCur = regions.getLabelRegion(node).localizingCursor();

                        // Für die Berechnung des Schwerpunkts
                        double[] positionSum = new double[labeling.numDimensions()];
                        int[] pixelCount = new int[labeling.numDimensions()];
                        Arrays.fill(pixelCount, 0);

                        // Iteration über Knotenpixel
                        while (nodeCur.hasNext()) {
                                nodeCur.fwd();

                                // Inputposition auf labelpixel setzen
                                nodeCur.localize(nodePixelPosition);

                                // Schwerpunktsberechnung: Addition aller Pixelcoordinaten/Anzahl der Pixel
                                for (int i = 0; i < labeling.numDimensions(); i++) {
                                        positionSum[i] += nodePixelPosition[i];
                                        pixelCount[i]++;
                                }

                                // ROI um nodepixel legen                        
                                roiAccess.setPosition(nodeCur);

                                Cursor<BitType> roiCursor = roiAccess.get().cursor();

                                // über ROI um Knotenpixel iterieren
                                while (roiCursor.hasNext()) {
                                        roiCursor.fwd();

                                        // wenn der Pixel weiß ist = eine Linie ist ...
                                        if (roiCursor.get().getRealDouble() != 0) {
                                                roiCursor.localize(nextOnLine);

                                                // ... die aber nicht zur current node gehört ...
                                                labelingRndAccess.setPosition(nextOnLine);
                                                if (!labelingRndAccess.get().contains(node)) {

                                                        // ... dann laufe weiter, bis ein weiterer knoten entdeckt wurde, der nicht er selbst ist.
                                                        Integer nextNode = walk(nextOnLine, nodePixelPosition, labelingRndAccess,
                                                                        regions.getExistingLabels());
                                                        if ((nextNode != null) && nextNode != node) {
                                                                nodeArray[node].add(new IntCell(nextNode));
                                                                //System.out.println("Nachbar von " + node + " ist " + nextNode);
                                                        }

                                                }
                                        }
                                }
                        }

                        // Schwerpunktsberechnung: Addition aller Pixelcoordinaten/Anzahl der Pixel
                        for (int i = 0; i < labeling.numDimensions(); i++) {
                                // Center berechnen + speichern
                                centers[node][i] = positionSum[i] / pixelCount[i];
                                //System.out.println("Node " + node + "[" + i + "]: " + centers[node][i]);
                        }

                }

                System.out.println("File1: Nodepositions");
                // Knoteniteration. Das erste Label ist 1, nicht 0...
                for (int i = 1; i < centers.length; i++) {

                        // über alle Nachbarn Iterieren
                        for (int neigh = 0; neigh < nodeArray[i].size(); neigh++) {
                                String outStr = "";

                                // dimensionsiteration 
                                for (int dim = 0; dim < labeling.numDimensions(); dim++) {
                                        if (dim != labeling.numDimensions() - 1) {
                                                outStr += String.format((Locale) null, "%.6f", centers[i][dim]) + "\t";
                                        } else {
                                                outStr += String.format((Locale) null, "%.6f", centers[i][dim]);
                                        }
                                }
                                System.out.println(outStr); //"Node " + i + ": " + 
                        }
                }

                System.out.println("File2: Neighbourpositions");
                // Knoteniteration. Das erste Label ist 1, nicht 0...
                for (int i = 1; i < centers.length; i++) {

                        // über alle Nachbarn Iterieren
                        for (int neigh = 0; neigh < nodeArray[i].size(); neigh++) {
                                String outStr = "";

                                // dimensionsiteration 
                                for (int dim = 0; dim < labeling.numDimensions(); dim++) {
                                        if (dim != labeling.numDimensions() - 1) {
                                                outStr += String.format((Locale) null, "%.6f", centers[nodeArray[i].get(neigh).getIntValue()][dim])
                                                                + "\t";
                                        } else {
                                                outStr += String.format((Locale) null, "%.6f", centers[nodeArray[i].get(neigh).getIntValue()][dim]);
                                        }
                                }

                                System.out.println(outStr);
                        }

                        //System.out.println(i + "\t" + nodeArray[i].size() + "\t" + nodeArray[i]);

                }

                // TODO: Ausgabe in Tabelle - end 2015 (maybe)
                // until then: copy data from logfile
                // the DataTableSpec of the final table
                // From new node Wizard (https://tech.knime.org/execute-0):
                /*
                DataTableSpec spec = new DataTableSpec(new DataColumnSpecCreator("Node", IntCell.TYPE).createSpec(), new DataColumnSpecCreator(
                                "Number Next Neighbours", IntCell.TYPE).createSpec(), new DataColumnSpecCreator("Next Neighbours",
                                ListCell.getCollectionType(IntCell.TYPE)).createSpec());
                
                BufferedDataContainer buf = exec.createDataContainer(spec);
                
                DataCell[] cells = new DataCell[3];
                
                for (int i = 1; i < nodeArray.length; i++) {
                        cells[0] = new IntCell(i);
                        cells[1] = new IntCell(nodeArray[i].size());
                        cells[2] = CollectionCellFactory.createListCell(nodeArray[i]);
                        buf.addRowToTable(new DefaultRow(("Node " + i), cells));
                }
                buf.close();
                BufferedDataTable outTable = buf.getTable();
                 */
        }

        // Läuft solange auf einer Linie weiter, bis der nächste Schritt gelabelt ist, gibt diesen Knoten zurück oder null wenn es eine Sackgasse ist
        private Integer walk(int[] current, int[] last, RandomAccess<LabelingType<Integer>> labelingRndAccess, Set<Integer> nodes) {

                // den nächsten Schritt ermitteln
                int[] next = new int[labeling.numDimensions()];

                // 3^n ROI generieren und um current legen TODO: roi & roiaccess als private variable
                final RectangleShape roi = new RectangleShape(1, false);
                RandomAccess<Neighborhood<BitType>> roiAccess = roi.neighborhoodsRandomAccessible(Views.extendValue(input, new BitType(false)))
                                .randomAccess();

                roiAccess.setPosition(current);
                Cursor<BitType> roiCursor = roiAccess.get().cursor();

                // über ROI iterieren
                while (roiCursor.hasNext()) {
                        roiCursor.fwd();

                        // Wenn der Pixel weiß ist...
                        if (roiCursor.get().getRealDouble() != 0) {
                                roiCursor.localize(next);

                                // Abbruchbedingungen checken

                                // Abbruchbedingung: wenn der nächste Schritt nicht zurückgeht...
                                if (!Arrays.equals(next, last) && !Arrays.equals(next, current)) {
                                        // ...über alle Knoten iterieren, und wenn next gelabelt ist diesen Knoten zurückgeben
                                        for (Integer node : nodes) {
                                                labelingRndAccess.setPosition(next);
                                                if (labelingRndAccess.get().contains(node)) {
                                                        return node;
                                                }
                                        }
                                        // wenn ein weißer Pixel gefunden ist, der weiß ist und nicht gelabelt, weitergehen
                                        return walk(next, current, labelingRndAccess, nodes);
                                }
                        }
                }
                // Abbruchbedingung: kein Pixel um current außer last ist mehr weiß: kein Nachbar
                return null;
        }
}
