import numpy as np
from Delaunator import Delaunator
import time

start_time = time.process_time()

#create random 2000 2d points between 4000000 and 4500000
points = np.random.randint(4000000,4500000,size=(1000,2))

point_creation_time = time.process_time() - start_time

triangles = Delaunator(points).triangles
halfedges = Delaunator(points).halfedges
hull = Delaunator(points).hull
coords = Delaunator(points).coords

print("\ncoords = "+str(coords))
print("\ntriangles = "+str(triangles))
print("\nhalfedges = "+str(halfedges))
print("\nhull = "+str(hull))

triangulation_time = time.process_time() - point_creation_time
print("\n" + str(point_creation_time))
print("\n" + str(triangulation_time))
