from Delaunator import *

points = []

for x in range (40,60,5):
    for y in range (11,20,3):
        points.append([x,y])

print(points)
triangles = Delaunator().run(points)
print(triangles)
