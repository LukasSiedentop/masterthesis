package org.knime.knip.pointextraction;

import java.util.Arrays;
import java.util.Iterator;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.labeling.Labeling;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

@SuppressWarnings("deprecation")
@Plugin(menu = {@Menu(label = "PointExtraction"), @Menu(label = "Nearest Neighbours")}, description = "Gets the Neighbours of Nodes which are connected by a pixel line. Takes a skeleton BitType Image and a Labeling with the Nodes.", headless = true, type = Command.class)
public class GetNearestNeighbours<BitType extends RealType<BitType>> implements Command {

        @Parameter(type = ItemIO.INPUT)
        private ImgPlus<BitType> input;

        @Parameter(type = ItemIO.INPUT)
        private Labeling<Integer> labeling;

        @Parameter(type = ItemIO.OUTPUT)
        private ImgPlus<UnsignedByteType> output;

        @Override
        public void run() {

                // javadoc.imagej.net/ImgLib2/net/imglib2/labeling/Labeling.html
                Iterator<Integer> nodes = labeling.getLabels().iterator();

                RandomAccess<BitType> inRndAccess = input.randomAccess();

                int[] nodePixelPosition = new int[labeling.numDimensions()];
                int[] nextOnLine = new int[labeling.numDimensions()];

                Integer node = 0;

                // ROI generieren mit Fenstergröße um den Ursprung 0^d
                double[] extend = new double[labeling.numDimensions()];
                Arrays.fill(extend, 3);
                double[] displacement = new double[labeling.numDimensions()];
                Arrays.fill(displacement, -1);
                final RectangleRegionOfInterest roi = new RectangleRegionOfInterest(new double[labeling.numDimensions()], extend);

                // Iteration über alle Knoten
                while (nodes.hasNext()) {
                        node = nodes.next();

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
                                                        System.out.println("Nachbar von " + node + " ist " + nextNode);
                                                }
                                        }
                                }
                                /*
                                                                System.out.println("Labelpixel: " + nodeCur.get());
                                                                System.out.println("Inputpixel: " + inRndAccess.get());
                                                                System.out.println("Position: " + Arrays.toString(posi));
                                                                System.out.println("Area: " + labeling.getArea(node));

                                                                RandomAccess<LabelingType<Integer>> labelingCur = labeling.randomAccess(input);
                                                                labelingCur.setPosition(posi);

                                                                System.out.println("LabelingCur: " + labelingCur.get());
                                                                */
                        }

                }

                System.out.println("feadsch");

                /* From new node Wizard (https://tech.knime.org/execute-0):
                int nrRows = ...; 
                int nrColumns = ...;
                BufferedDataContainer buf = exec.createBufferedDataContainer(spec);
                for (int j = 0; j < nrRows; j++) {
                DataCell[] cells = new DataCell[nrColumns];
                for (int i = 0; i < nrColumns; i++) {
                cells[i] = new DoubleCell(i * Math.PI);
                }
                DataRow row = new DefaultRow(
                new StringCell(“RowKey_” + j, cells);
                buf.addRowToTable(row);
                }
                buf.close();
                BufferedDataTable table = buf.getTable();
                 */

        }

        private Integer walk(int[] current, int[] last) {

                // den nächsten Schritt ermitteln
                int[] next = new int[labeling.numDimensions()];

                // ROI generieren mit Fenstergröße um den Ursprung 0^d
                double[] extend = new double[labeling.numDimensions()];
                Arrays.fill(extend, 3);
                double[] displacement = new double[labeling.numDimensions()];
                Arrays.fill(displacement, -1);
                final RectangleRegionOfInterest roi = new RectangleRegionOfInterest(new double[labeling.numDimensions()], extend);
                // ROI um current legen
                roi.setOrigin(copyFromIntArray(current));
                roi.move(displacement);

                Cursor<BitType> roiCursor = roi.getIterableIntervalOverROI(input).cursor();

                // über ROI iterieren
                while (roiCursor.hasNext()) {
                        roiCursor.fwd();

                        if (roiCursor.get().getRealDouble() != 0) {
                                roiCursor.localize(next);

                                // Abbruchbedingungen checken
                                // Abbruchbedingung: wenn der nächste Schritt nicht zurückgeht...
                                if (!Arrays.equals(next, last) && !Arrays.equals(next, current)) {
                                        Iterator<Integer> nodesIter = labeling.getLabels().iterator();
                                        Integer node = 0;
                                        // über alle Knoten iterieren, und wenn next gelabelt ist diesen Knoten zurückgeben
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
                return -99;

        }

        // Hack: convert to int[] double[], damit roi.contains funktioniert. Diskussion:
        // http://stackoverflow.com/questions/12729139/copy-contents-of-an-int-array-to-a-double-array-in-java
        private static double[] copyFromIntArray(int[] source) {
                double[] dest = new double[source.length];
                for (int i = 0; i < source.length; i++) {
                        dest[i] = source[i];
                }
                return dest;
        }
}
