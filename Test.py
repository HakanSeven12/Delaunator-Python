import numpy as np
from Delaunator import Delaunator

#create random 2000 2d points between 4000000 and 4500000
points = np.random.randint(4000000,4500000,size=(2000,2))

triangles = Delaunator().run(points)
print(triangles)
