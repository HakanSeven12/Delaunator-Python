# Delaunator-Python
Fast Delaunay triangulation of 2D points implemented in Python.

This code was ported from [Mapbox's Delaunator Project](https://github.com/mapbox/delaunator) (JavaScript).

# Usage
```python
from Delaunator import Delaunator

points = [[168, 180], [168, 178], [168, 179], [168, 181], [168, 183], ...]

triangles = Delaunator().run(points)
print(triangles)
```
