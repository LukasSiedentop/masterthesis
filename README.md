# Masterthesis

Stuff for my masterthesis.

## KNIP Nodes
[Javadoc](http://javadoc.imagej.net/)

### Neighbours Count
Counts the direct neighbours of each white Pixel in a BitType Image (assuming 8-connectivity). Returns the counting value as a ByteType Image, excluding the white Pixel itself. Useful to find intersections in a skeleton BitType Image.

### Connected Neighbours
Given a (skeleton) BitType Image and a Labeling labeling the intersections, this node evaluates the connected neighbours of each intersection. An intersection is neighbouring another if it's "connected" via a white pixel trail. It prints the positions of the Nodes and Neighbours on the Console. BufferedTable support to the ImageJ wrapper will (maybe) be implemented end 2015.

### Shift/Align/Translate/Shear
Shears an image. Each pixel is shifted in a given direction with a distance depending on the position in a chosen dimension and a given factor.

## Shapes
Little python script that generates `.gwl` files containing coordinates representing a structure (e.g. a woodpile, trampoline, ...) which can be fed to the NanoScribe unit.

## Spectrometer
Little python script that enables a measurement on our optical table. Basically, it controls a monochromator (SpectraPro-2150i), a stage motor (ESP100) and reads out the amplified signal of a photodiode via a data acquisition device (NI USB-6211).

## Point Pattern Tool
C++ tool that generates various statistics of point patterns, such as the distributions of lengths or angles. Uses the [gnuplot-iostream interface](http://www.stahlke.org/dan/gnuplot-iostream/)(source included) and includes [voro++](http://math.lbl.gov/voro++/) as a static library.
