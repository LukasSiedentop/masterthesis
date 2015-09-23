package org.knime.knip.pointextraction;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import net.imagej.ImgPlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.NativeType;

@SuppressWarnings("deprecation")
@Plugin(menu = {@Menu(label = "PointExtraction"),
                @Menu(label = "Shift/Align/Translate/Shear")}, description = "Shifts the Image in a chosen direction as a linear function of the second chose dimension (e.g. shift in Y-direction with step * Z-position). Attention: no interpolation is taking place!", headless = true, type = Command.class)
public class ShiftYwithZ<T extends NativeType<T>> implements Command {

        @Parameter(type = ItemIO.INPUT)
        private ImgPlus<T> input;

        @Parameter(type = ItemIO.INPUT, label = "step", description = "Distance to shift: step * position in chosen dimension")
        private double step = 1;

        // TODO: dimension selection
        @Parameter(type = ItemIO.INPUT, label = "shifting direction", description = "Direction in which is to be shifted.")
        private int shiftDimension = 1;

        // TODO: dimension selection
        @Parameter(type = ItemIO.INPUT, label = "shifting distance", description = "Distance to shift with: step * position in this dimension.")
        private int argumentDimension = 2;

        @Parameter(type = ItemIO.OUTPUT)
        private ImgPlus<T> output;

        @Override
        public void run() {
                // create empty output image with adjusted size
                long[] dims = new long[input.numDimensions()];
                input.dimensions(dims);
                dims[shiftDimension] += dims[argumentDimension] * step;
                output = new ImgPlus<T>(new ArrayImgFactory<T>().create(dims, input.firstElement()));

                // Hack to get dimension names correctly
                for (int d = 0; d < input.numDimensions(); d++) {
                        output.setAxis(input.axis(d), d);
                }

                // access the pixels of the output image
                final RandomAccess<T> outAccess = output.randomAccess();

                // cursor over input image
                final Cursor<T> inCursor = input.localizingCursor();

                while (inCursor.hasNext()) {
                        inCursor.fwd();

                        // set outaccess on position of incursor plus shift
                        outAccess.setPosition(inCursor);
                        outAccess.move((int) (inCursor.getIntPosition(argumentDimension) * step), shiftDimension);

                        // set output pixel
                        outAccess.get().set(inCursor.get());
                }
        }
}
