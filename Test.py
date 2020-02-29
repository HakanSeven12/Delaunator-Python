import numpy as np
from Delaunator import Delaunator
import time

start_time = time.clock()

#create random 2000 2d points between 4000000 and 4500000
points = np.random.randint(4000000,4500000,size=(100000,2))

point_creation_time = time.clock() - start_time

triangles = Delaunator().run(points)
print(triangles)

triangulation_time = time.clock() - point_creation_time
print(point_creation_time)
print(triangulation_time)
