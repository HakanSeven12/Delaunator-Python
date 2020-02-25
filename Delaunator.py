EPSILON = 2**-52
EDGE_STACK =[]

import math

class Delaunator:
    def constructor(self, coords):
        if (len(points) < 3):
            raise ValueError("Need at least 3 points")
        self.coords = coords
        n = len(self.coords) >> 1

        # arrays that will store the triangulation graph
        maxTriangles = max(2 * n - 5, 0)
        self._triangles = [None] * maxTriangles * 3
        self._halfedges = [None] * maxTriangles * 3

        # temporary arrays for tracking the edges of the advancing convex hull
        self.hashSize = math.ceil(math.sqrt(n))
        self.hullPrev = [None] * n # edge to prev edge
        self.hullNext = [None] * n # edge to next edge
        self.hullTri = [None] * n # edge to adjacent triangle
        self.hullHash = [-1] * self.hashSize # angular edge hash

        # temporary arrays for sorting points
        self._ids =  [None] * n
        self._dists = [None] * n

        self.update()

    def update(self):
        n = len(self.coords) >> 1

        # populate an array of point indices; calculate input data bbox
        minX = math.inf
        minY = math.inf
        maxX = -math.inf
        maxY = -math.inf

        for i in range(0,n):
            x = coords[2 * i]
            y = coords[2 * i + 1]
            if (x < minX): minX = x
            if (y < minY): minY = y
            if (x > maxX): maxX = x
            if (y > maxY): maxY = y
            self._ids[i] = i

        cx = (minX + maxX) / 2
        cy = (minY + maxY) / 2

        minDist = math.inf
        i0 = 0
        i1 = 0
        i2 = 0

        # pick a seed point close to the center
        for i in range(0,n):
            d = dist(cx, cy, coords[2 * i], coords[2 * i + 1])
            if (d < minDist):
                i0 = i
                minDist = d

        i0x = coords[2 * i0]
        i0y = coords[2 * i0 + 1]

        minDist = math.inf

        # find the point closest to the seed
        for i in range(0,n):
            if (i == i0):
                continue

            d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1])
            if (d < minDist and d > 0):
                i1 = i
                minDist = d

        i1x = coords[2 * i1]
        i1y = coords[2 * i1 + 1]

        minRadius = math.inf

        # find the third point which forms the smallest circumcircle with the first two
        for i in range(0,n):
            if (i == i0 or i == i1):
                continue

            r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1])
            if (r < minRadius):
                i2 = i
                minRadius = r

        i2x = coords[2 * i2]
        i2y = coords[2 * i2 + 1]

        if (minRadius == math.inf):
            # order collinear points by dx (or dy if all x are identical)
            # and return the list as a hull
            for i in range(0,n):
                self._dists[i] = (coords[2 * i] - coords[0]) or (coords[2 * i + 1] - coords[1])

            quicksort(self._ids, self._dists, 0, n - 1)
            hull =  (n)
            j = 0
            d0 = -math.inf
            for i in range(0,n):
                id = self._ids[i]
                if (self._dists[id] > d0):
                    j=j+1
                    hull[j] = id
                    d0 = self._dists[id]

            self.hull = hull.subarray(0, j)
            self.triangles =  (0)
            self.halfedges =  (0)
            return

        # swap the order of the seed points for counter-clockwise orientation
        if (orient(i0x, i0y, i1x, i1y, i2x, i2y)):
            i = i1
            x = i1x
            y = i1y
            i1 = i2
            i1x = i2x
            i1y = i2y
            i2 = i
            i2x = x
            i2y = y

        center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y)
        self._cx = center[0]
        self._cy = center[1]

        for i in range(0,n):
            self._dists[i] = dist(coords[2 * i], coords[2 * i + 1], center[0], center[1])

        # sort the points by distance from the seed triangle circumcenter
        quicksort(self._ids, self._dists, 0, n - 1)

        # set up the seed triangle as the starting hull
        self._hullStart = i0
        hullSize = 3

        self.hullNext[i0] = self.hullPrev[i2] = i1
        self.hullNext[i1] = self.hullPrev[i0] = i2
        self.hullNext[i2] = self.hullPrev[i1] = i0

        self.hullTri[i0] = 0
        self.hullTri[i1] = 1
        self.hullTri[i2] = 2

        self.hullHash[self._hashKey(i0x, i0y)] = i0
        self.hullHash[self._hashKey(i1x, i1y)] = i1
        self.hullHash[self._hashKey(i2x, i2y)] = i2

        self.trianglesLen = 0
        self._addTriangle(i0, i1, i2, -1, -1, -1)
        xp=0
        yp=0

        for k in range(0,len(self._ids)):
            i = self._ids[k]
            x = coords[2 * i]
            y = coords[2 * i + 1]

            # skip near-duplicate points
            if (k > 0 and math.abs(x - xp) <= EPSILON and math.abs(y - yp) <= EPSILON):
                xp = x
                yp = y

            # skip seed triangle points
            if (i == i0 or i == i1 or i == i2):
                # find a visible edge on the convex hull using edge hash
                start = 0
            key = self._hashKey(x, y)

            for j in range(0,self.hashSize):
                start = self.hullHash[(key + j) % self.hashSize]
                if (start != -1 and start != self.hullNext[start]):
                    break

            start = self.hullPrev[start]
            e = start
            q = self.hullNext[e]
            
            while True:
                if orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]):
                    break
                e = q
                if (e == start):
                    e = -1
                    break

            if (e == -1): # likely a near-duplicate point; skip it
                continue

            # add the first triangle from the point
            t = self._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e])

            # recursively flip triangles from the point until they satisfy the Delaunay condition
            hullTri[i] = self._legalize(t + 2)
            hullTri[e] = t # keep track of boundary triangles on the hull
            ++hullSize

            # walk forward through the hull, adding more triangles and flipping recursively
            n = hullNext[e]
            q = hullNext[n]
            
            while (orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])):
                t = self._addTriangle(n, i, q, hullTri[i], -1, hullTri[n])
                hullTri[i] = self._legalize(t + 2)
                hullNext[n] = n # mark as removed
                hullSize-=1
                n = q

            # walk backward from the other side, adding more triangles and flipping
            if (e == start):
                q = hullPrev[e]
                while (orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])):
                    t = self._addTriangle(q, i, e, -1, hullTri[e], hullTri[q])
                    self._legalize(t + 2)
                    hullTri[q] = t
                    hullNext[e] = e # mark as removed
                    hullSize-=1
                    e = q

            # update the hull indices
            self._hullStart = hullPrev[i] = e
            hullNext[e] = hullPrev[n] = i
            hullNext[i] = n

            # save the two new edges in the hash table
            hullHash[self._hashKey(x, y)] = i
            hullHash[self._hashKey(coords[2 * e], coords[2 * e + 1])] = e

        self.hull = (hullSize)
        e = self._hullStart
        for i in range(0,hullSize):
            self.hull[i] = e
            e = hullNext[e]

        # trim typed triangle mesh arrays
        self.triangles = self._triangles.subarray(0, self.trianglesLen)
        self.halfedges = self._halfedges.subarray(0, self.trianglesLen)
        print (self.triangles)

    def _hashKey(self,x, y):
        return math.floor(pseudoAngle(x - self._cx, y - self._cy) * self.hashSize) % self.hashSize

    def _legalize(a):
        i = 0
        ar = 0

        # recursion eliminated with a fixed-size stack
        while True:
            b = halfedges[a]
            """
              if the pair of triangles doesn't satisfy the Delaunay condition
              (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
              then do the same check/flip recursively for the new pair of triangles
             
                        pl                    pl
                       /||\                  /  \
                    al/ || \bl            al/    \a
                     /  ||  \              /      \
                    /  a||b  \    flip    /___ar___\
                  p0\   ||   /p1   =>   p0\---bl---/p1
                     \  ||  /              \      /
                    ar\ || /br             b\    /br
                       \||/                  \  /
                        pr                    pr
             
            """
            a0 = a - a % 3
            ar = a0 + (a + 2) % 3

            if (b == -1): # convex hull edge
                if (i == 0):
                    break
                a = EDGE_STACK[--i]
                continue

            b0 = b - b % 3
            al = a0 + (a + 1) % 3
            bl = b0 + (b + 2) % 3

            p0 = triangles[ar]
            pr = triangles[a]
            pl = triangles[al]
            p1 = triangles[bl]

            illegal = inCircle(
                coords[2 * p0], coords[2 * p0 + 1],
                coords[2 * pr], coords[2 * pr + 1],
                coords[2 * pl], coords[2 * pl + 1],
                coords[2 * p1], coords[2 * p1 + 1])

            if (illegal):
                triangles[a] = p1
                triangles[b] = p0

                hbl = halfedges[bl]

                # edge swapped on the other side of the hull (rare); fix the halfedge reference
                if (hbl == -1):
                    e = self._hullStart
                    if (self._hullTri[e] == bl):
                        self._hullTri[e] = a
                        break

                    e = self._hullPrev[e]
                    while (e != self._hullStart):
                        if (self._hullTri[e] == bl):
                            self._hullTri[e] = a
                            break

                self._link(a, hbl)
                self._link(b, halfedges[ar])
                self._link(ar, bl)

                br = b0 + (b + 1) % 3

                # don't worry about hitting the cap: it can only happen on extremely degenerate input
                if (i < EDGE_STACK.length):
                    EDGE_STACK[++i] = br

            else:
                if (i == 0):
                   break
                a = EDGE_STACK[--i]

        return ar

    def _link(self,a, b):
        self._halfedges[a] = b
        if (b != -1):
            self._halfedges[b] = a

    # add a new triangle given vertex indices and adjacent half-edge ids
    def _addTriangle(self,i0, i1, i2, a, b, c):
        t = self.trianglesLen

        self._triangles[t] = i0
        self._triangles[t + 1] = i1
        self._triangles[t + 2] = i2

        self._link(t, a)
        self._link(t + 1, b)
        self._link(t + 2, c)

        self.trianglesLen += 3

        return t

# monotonically increases with real angle, but doesn't need expensive trigonometry
def pseudoAngle(dx, dy):
    try:
      p = dx / (abs(dx) + abs(dy))
    except:
      p = 0

    if (dy > 0):
        return (3 - p) / 4 # [0..1]
    else:
        return (1 + p) / 4 # [0..1]

def dist(ax, ay, bx, by):
    dx = ax - bx
    dy = ay - by
    return dx * dx + dy * dy

# return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
def orientIfSure(px, py, rx, ry, qx, qy):
    l = (ry - py) * (qx - px)
    r = (rx - px) * (qy - py)
    
    if (abs(l - r) >= 3.3306690738754716e-16 * abs(l + r)):
        return l - r
    else:
        return 0

# a more robust orientation test that's stable in a given triangle (to fix robustness issues)
def orient(rx, ry, qx, qy, px, py):
    return (orientIfSure(px, py, rx, ry, qx, qy) or
        orientIfSure(rx, ry, qx, qy, px, py) or
        orientIfSure(qx, qy, px, py, rx, ry)) < 0

def inCircle(ax, ay, bx, by, cx, cy, px, py):
    dx = ax - px
    dy = ay - py
    ex = bx - px
    ey = by - py
    fx = cx - px
    fy = cy - py

    ap = dx * dx + dy * dy
    bp = ex * ex + ey * ey
    cp = fx * fx + fy * fy

    return dx * (ey * cp - bp * fy) -\
           dy * (ex * cp - bp * fx) +\
           ap * (ex * fy - ey * fx) < 0

def circumradius(ax, ay, bx, by, cx, cy):
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    try:
      d = 0.5 / (dx * ey - dy * ex)
    except:
      d = 0

    x = (ey * bl - dy * cl) * d
    y = (dx * cl - ex * bl) * d

    return x * x + y * y

def circumcenter(ax, ay, bx, by, cx, cy):
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    try:
        d = 0.5 / (dx * ey - dy * ex)
    except:
        d = 0

    x = ax + (ey * bl - dy * cl) * d
    y = ay + (dx * cl - ex * bl) * d

    return x, y

def quicksort(ids, dists, left, right):
    if (right - left <= 20):
        for i in range(left + 1,right):
            temp = ids[i]
            tempDist = dists[temp]
            j = --i
            while (j >= left and dists[ids[j]] > tempDist):
                ids[j + 1] = ids[--j]
            ids[j + 1] = temp;

    else:
        median = (left + right) >> 1
        i = left + 1
        j = right
        swap(ids, median, i)

        if (dists[ids[left]] > dists[ids[right]]):
            swap(ids, left, right)

        if (dists[ids[i]] > dists[ids[right]]):
            swap(ids, i, right)

        if (dists[ids[left]] > dists[ids[i]]):
            swap(ids, left, i)

        temp = ids[i]
        tempDist = dists[temp]

        while True:
            ++i
            while (dists[ids[i]] < tempDist):
                ++i
            --j
            while (dists[ids[j]] > tempDist):
                --j

            if (j < i):
                break

            swap(ids, i, j);

        ids[left + 1] = ids[j];
        ids[j] = temp;

        if (right - i + 1 >= j - left):
            quicksort(ids, dists, i, right)
            quicksort(ids, dists, left, j - 1)

        else:
            quicksort(ids, dists, left, j - 1)
            quicksort(ids, dists, i, right)

def swap(arr, i, j):
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp

points = [[168, 180], [168, 178], [168, 179], [168, 181], [168, 183]]
n = len(points)
coords = [None] * n * 2

for i in range(0,n):
    p = points[i]
    coords[2 * i] = (p[0])
    coords[2 * i+1] = (p[1])

Delaunator().constructor(coords)
