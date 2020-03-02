# Delaunator-Python
Fast Delaunay triangulation of 2D points implemented in Python.

This code was ported from [Mapbox's Delaunator Project](https://github.com/mapbox/delaunator) (JavaScript).

Test on FreeCAD.

![Test](https://user-images.githubusercontent.com/3831435/75558770-60245280-5a53-11ea-8d1f-855c9e3f9c13.png)

## Usage
```python
from Delaunator import Delaunator

points = [[168, 180], [168, 178], [168, 179], [168, 181], [168, 183], ...]

triangles = Delaunator(points).triangles
print(triangles)

>> [623, 636, 619,  636, 444, 619, ...]
```

## API Reference

### Delaunator(points)

Constructs a delaunay triangulation object given an array of points (`[x, y]` by default).
Duplicate points are skipped.

### Delaunator(points).triangles

An array of triangle vertex indices (each group of three numbers forms a triangle).
All triangles are directed counterclockwise.

To get the coordinates of all triangles, use:

```python
coordinates = []

for i in range(0, len(triangles), 3):
    coordinates.append([
        points[triangles[i]],
        points[triangles[i + 1]],
        points[triangles[i + 2]]])
```

### Delaunator(points).halfedges

An array of triangle half-edge indices that allows you to traverse the triangulation.
`i`-th half-edge in the array corresponds to vertex `triangles[i]` the half-edge is coming from.
`halfedges[i]` is the index of a twin half-edge in an adjacent triangle
(or `-1` for outer half-edges on the convex hull).

The flat array-based data structures might be counterintuitive,
but they're one of the key reasons this library is fast.

### Delaunator(points).hull

An array of indices that reference points on the convex hull of the input data, counter-clockwise.

### Delaunator(points).coords

An array of input coordinates in the form `[x0, y0, x1, y1, ....]`,
of the type provided in the constructor.

### Delaunator(points).update()

Updates the triangulation if you modified `Delaunator(points).coords` values in place, avoiding expensive memory allocations.
Useful for iterative relaxation algorithms such as [Lloyd's](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm).
