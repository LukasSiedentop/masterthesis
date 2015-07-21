# Masterthesis

Stuff for my masterthesis.

## KNIP Nodes
[Javadoc](http://javadoc.imagej.net/)

### Neighbours Count
Counts the direct neighbours of each white Pixel in a BitType Image (assuming 8-connectivity). Returns the counting value as a ByteType Image, excluding the white Pixel itself. Useful to find intersections in a skeleton BitType Image.

### Connected Neighbours
Given a (skeleton) BitType Image and a Labeling labeling the intersections, this node evaluates the connected neighbours of each intersection. An intersection is neighbouring another if it's "connected" via a white pixel trail.

Example (1-6 = nodes; - / \ = connections):

```
       5---6
      / \
     /   \
    /     3--4
1--2
```

Would result in a table:

```
Node	Number of Connections	Connected Nodes
1		1		\[2\]
2		2		\[1, 5\]
3		2		\[5, 4\]
4		1		\[3\]
5		3		\[2, 3, 6\]
6		1		\[5\]
```

## Shapes
Little python script that generates `.gwl` files containing coordinates representing a structure (e.g. a woodpile, trampoline, ...) which can be fed to the NanoScribe unit.


