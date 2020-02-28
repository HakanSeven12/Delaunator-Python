# Delaunator-Python
Fast Delaunay triangulation of 2D points implemented in Python.

This code was ported from [Mapbox's Delaunator Project](https://github.com/mapbox/delaunator) (JavaScript).

Test on FreeCAD.

![Test](https://user-images.githubusercontent.com/3831435/75558770-60245280-5a53-11ea-8d1f-855c9e3f9c13.png)

# Usage
```python
from Delaunator import Delaunator

points = [[168, 180], [168, 178], [168, 179], [168, 181], [168, 183], ...]

triangles = Delaunator().run(points)
print(triangles)

>> [623, 636, 619,  636, 444, 619, ...]
```
