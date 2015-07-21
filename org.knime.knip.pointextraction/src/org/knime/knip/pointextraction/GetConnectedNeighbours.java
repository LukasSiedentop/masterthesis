package org.knime.knip.pointextraction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

import org.knime.core.data.def.IntCell;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.labeling.Labeling;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;

@SuppressWarnings("deprecation")
@Plugin(menu = {@Menu(label = "PointExtraction"),
                @Menu(label = "Connected Neighbours")}, description = "Gets the Neighbours of Nodes which are connected by a pixel line. Takes a skeleton BitType Image and a Labeling with the Nodes.", headless = true, type = Command.class)
public class GetConnectedNeighbours<BitType extends RealType<BitType>> implements Command {

        @Parameter(type = ItemIO.INPUT)
        private ImgPlus<BitType> input;

        @Parameter(type = ItemIO.INPUT)
        private Labeling<Integer> labeling;

        // geht nicht
        //@Parameter(type = ItemIO.INPUT)
        //private ExecutionContext exec;

        //@Parameter(type = ItemIO.OUTPUT)
        //private BufferedDataTable output;

        @Parameter(type = ItemIO.OUTPUT)
        private ImgPlus<UnsignedByteType> output;

        @Override
        public void run() {
                RandomAccess<BitType> inRndAccess = input.randomAccess();

                Iterator<Integer> nodes = labeling.getLabels().iterator();
                Integer node = 0;

                // 3^n ROI generieren
                final RectangleRegionOfInterest roi = getROI();
                double[] displacement = new double[labeling.numDimensions()];
                Arrays.fill(displacement, -1);

                // Positionsarrays
                int[] nodePixelPosition = new int[labeling.numDimensions()];
                int[] nextOnLine = new int[labeling.numDimensions()];

                // Datenstruktur um Ergebnis aufzunehmen
                @SuppressWarnings("unchecked")
                ArrayList<IntCell>[] nodeArray = ((ArrayList<IntCell>[]) new ArrayList[labeling.getLabels().size() + 1]);

                // Iteration über alle Knoten
                while (nodes.hasNext()) {
                        node = nodes.next();

                        nodeArray[node] = new ArrayList<IntCell>(0);

                        // Cursor über Knotenpixel
                        Cursor<BitType> nodeCur = labeling.getIterableRegionOfInterest(node).getIterableIntervalOverROI(input).localizingCursor();

                        // Iteration über Knotenpixel
                        while (nodeCur.hasNext()) {
                                nodeCur.fwd();

                                // Inputposition auf labelpixel setzen
                                nodeCur.localize(nodePixelPosition);
                                inRndAccess.setPosition(nodePixelPosition);

                                // ROI um nodepixel legen
                                roi.setOrigin(nodeCur);
                                roi.move(displacement);

                                Cursor<BitType> roiCursor = roi.getIterableIntervalOverROI(input).cursor();

                                // über ROI um Knotenpixel iterieren
                                while (roiCursor.hasNext()) {
                                        roiCursor.fwd();

                                        if (roiCursor.get().getRealDouble() != 0) {
                                                roiCursor.localize(nextOnLine);

                                                if (!labeling.getRegionOfInterest(node).contains(copyFromIntArray(nextOnLine))) {

                                                        Integer nextNode = walk(nextOnLine, nodePixelPosition);
                                                        if (nextNode != null) {
                                                                nodeArray[node].add(new IntCell(nextNode));
                                                                //System.out.println("Nachbar von " + node + " ist " + nextNode);
                                                        }

                                                }
                                        }
                                }
                        }

                }

                System.out.println("Node\tNumber of Connections\tConnected Nodes");
                for (int i = 1; i < nodeArray.length; i++) {
                        System.out.println(i + "\t" + nodeArray[i].size() + "\t" + nodeArray[i]);
                }

                // TODO: Ausgabe in Tabelle: ExecutionContext - woher? Wie die outTable rausgeben?
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

        // Läuft solange auf einer Linie weiter, bis der nächste Schritt gelabelt ist, gibt diesen Knoten zurück oder -99 wenn es eine Sackgasse ist
        private Integer walk(int[] current, int[] last) {

                // den nächsten Schritt ermitteln
                int[] next = new int[labeling.numDimensions()];

                // 3^n ROI generieren und um current legen
                final RectangleRegionOfInterest roi = getROI();
                double[] displacement = new double[labeling.numDimensions()];
                Arrays.fill(displacement, -1);
                roi.setOrigin(copyFromIntArray(current));
                roi.move(displacement);

                Cursor<BitType> roiCursor = roi.getIterableIntervalOverROI(input).cursor();

                // über ROI iterieren
                while (roiCursor.hasNext()) {
                        roiCursor.fwd();

                        // Wenn der Pixel weiß ist...
                        if (roiCursor.get().getRealDouble() != 0) {
                                roiCursor.localize(next);

                                // Abbruchbedingungen checken

                                // Abbruchbedingung: wenn der nächste Schritt nicht zurückgeht...
                                if (!Arrays.equals(next, last) && !Arrays.equals(next, current)) {
                                        Iterator<Integer> nodesIter = labeling.getLabels().iterator();
                                        Integer node = 0;
                                        // ...über alle Knoten iterieren, und wenn next gelabelt ist diesen Knoten zurückgeben
                                        while (nodesIter.hasNext()) {
                                                node = nodesIter.next();
                                                if (labeling.getRegionOfInterest(node).contains(copyFromIntArray(next))) {
                                                        return node;
                                                }
                                        }

                                        // wenn ein weißer Pixel gefunden ist, der weiß ist und nicht gelabelt, weitergehen
                                        return walk(next, current);
                                }

                        }
                }

                // Abbruchbedingung: kein Pixel um current außer last ist mehr weiß: kein Nachbar
                return null;

        }

        // Hack: convert to int[] double[], damit roi.contains und roi.setOrigin funktioniert. Diskussion:
        // http://stackoverflow.com/questions/12729139/copy-contents-of-an-int-array-to-a-double-array-in-java
        private static double[] copyFromIntArray(int[] source) {
                double[] dest = new double[source.length];
                for (int i = 0; i < source.length; i++) {
                        dest[i] = source[i];
                }
                return dest;
        }

        // returns a rectangle roi of extend 3^d 
        private final RectangleRegionOfInterest getROI() {
                double[] extend = new double[labeling.numDimensions()];
                Arrays.fill(extend, 3);
                return new RectangleRegionOfInterest(new double[labeling.numDimensions()], extend);
        }
}
