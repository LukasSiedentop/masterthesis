# Masterthesis

Stuff for my masterthesis.

## KNIP Nodes
[Javadoc](http://javadoc.imagej.net/)

### Neighbours Count
Counts the direct neighbours of each white Pixel in a BitType Image (assuming 8-connectivity). Returns the counting value as a ByteType Image, excluding the white Pixel itself. Useful to find intersections in a skeleton BitType Image.

### Nearest Neighbours
Given a (skeleton) BitType Image and a Labeling labeling the intersections, this node evaluates the neighbours of each intersection. An intersection is neighbouring another if it's "connected" via a white pixel trail.

## Shapes
Little python script that generates `.gwl` files containing coordinates representing a structure (e.g. a woodpile, trampoline, ...) which can be fed to the NanoScribe unit.


