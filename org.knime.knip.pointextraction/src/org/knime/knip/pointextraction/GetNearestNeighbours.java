package org.knime.knip.pointextraction;

import java.util.Arrays;
import java.util.Iterator;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.labeling.Labeling;
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

                long[] posi = new long[labeling.numDimensions()];
                Integer node = 0;

                // Iteration über alle Knoten
                while (nodes.hasNext()) {
                        node = nodes.next();

                        // Cursor über Knotenpixel
                        Cursor<BitType> nodeCur = labeling.getIterableRegionOfInterest(node).getIterableIntervalOverROI(input).localizingCursor();

                        // Iteration über Knotenpixel
                        while (nodeCur.hasNext()) {
                                nodeCur.fwd();

                                // Inputposition auf labelpixel setzen
                                nodeCur.localize(posi);
                                inRndAccess.setPosition(posi);

                                System.out.println("Labelpixel: " + nodeCur.get());
                                System.out.println("Inputpixel: " + inRndAccess.get());
                                System.out.println("Position: " + Arrays.toString(posi));
                                System.out.println("Area: " + labeling.getArea(node));
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
}
