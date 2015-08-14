package org.knime.knip.pointextraction;

import java.util.Arrays;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.roi.RectangleRegionOfInterest;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;

@SuppressWarnings("deprecation")
@Plugin(menu = {@Menu(label = "PointExtraction"),
                @Menu(label = "Neighbours Count")}, description = "Counts the Neighbours of a white Pixel in a BitType Image and returns them in a ByteTypeImage.", headless = true, type = Command.class)
public class NeighboursCount<BitType extends RealType<BitType>> implements Command {

        @Parameter(type = ItemIO.INPUT)
        private ImgPlus<BitType> input;

        @Parameter(type = ItemIO.INPUT, label = "Window Size", description = "Windowsize in which neighbours are searched for.")
        private double windowSize = 3;

        @Parameter(type = ItemIO.OUTPUT)
        private ImgPlus<UnsignedByteType> output;

        @Override
        public void run() {
                // TODO: Knoten als Labeling ausgeben, wenn mal langweilig ist...

                // create empty output image (arrayimg of type bittype)
                output = new ImgPlus<UnsignedByteType>(new ArrayImgFactory<UnsignedByteType>().create(input, new UnsignedByteType()));

                // Hack um die Dimensionsbenamsung richtig zu setzen
                for (int d = 0; d < input.numDimensions(); d++) {
                        output.setAxis(input.axis(d), d);
                }

                // access the pixels of the output image
                final RandomAccess<UnsignedByteType> outAccess = output.randomAccess();

                // cursor over input image
                final Cursor<BitType> inCursor = input.localizingCursor();

                // ROI generieren mit Fenstergröße um den Ursprung 0^d
                double[] extend = new double[inCursor.numDimensions()];
                Arrays.fill(extend, windowSize);
                double[] displacement = new double[inCursor.numDimensions()];
                Arrays.fill(displacement, -windowSize / 2.0);
                final RectangleRegionOfInterest roi = new RectangleRegionOfInterest(new double[inCursor.numDimensions()], extend);

                while (inCursor.hasNext()) {
                        // weitergehen
                        inCursor.fwd();

                        // position setzen
                        int[] posi = new int[inCursor.numDimensions()];
                        inCursor.localize(posi);

                        // nicht tun, wenn Pixel selbst schwarz TODO: als option einbauen?
                        if (inCursor.get().getRealDouble() == 0) {
                                continue;
                        }

                        // set outaccess on position of incursor
                        outAccess.setPosition(inCursor);

                        // Fenster um den aktuellen Ort legen (Origin ist im linken oberen Eck)
                        roi.setOrigin(inCursor);
                        roi.move(displacement);

                        // Cursor im ROI bekommen
                        final Cursor<BitType> roiCursor = roi.getIterableIntervalOverROI(input).cursor();

                        // Nachbarn zusammenzählen und aufaddieren (sich selbst nicht mitzählen)
                        byte neighbours = -1;
                        while (roiCursor.hasNext()) {
                                roiCursor.fwd();

                                //outAccess.setPosition(roiCursor);
                                //outAccess.get().set((byte) 50);

                                if (roiCursor.get().getRealDouble() != 0)
                                        neighbours += 1;

                        }
                        outAccess.get().set(neighbours);
                }

        }
}
