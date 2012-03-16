import sys
import numpy
import math
import scipy.spatial.distance as distance

import triangle

with open('cat1.txt') as fd:
    cat1 = numpy.loadtxt(fd)

with open('cat2.txt') as fd:
    cat2 = numpy.loadtxt(fd)

def sort_by_mag(a):
    # http://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
    return a[a[:,2].argsort()]

cat1s = sort_by_mag(cat1)
cat2s = sort_by_mag(cat2)

def normalize1(a):
    mi = a.min()
    mx = a.max()

    a = (a - mi) / (mx - mi)
    return a

def normaliza(a):
    b = a[:]
    b[:,0] = normalize1(a[:,0])
    b[:,1] = normalize1(a[:,1])

    return b

reject_scale = 10.0

common = min(cat1s.shape[0], cat2s.shape[0])
print common

a = normaliza(cat1s[:common,:])
l = common
tl1 = []
for i in range(l):
    for j in range(i + 1, l):
        for k in range(j + 1, l):
            r = a[[i,j,k], :]
            tng = triangle.create_triang(r)
            # if scale R > reject_scale, reject
            if tng[5] < reject_scale:
                tl1.append(tng)

a = normaliza(cat2s[:common,:])
l = common
tl2 = []
for i in range(l):
    for j in range(i + 1, l):
        for k in range(j + 1, l):
            r = a[[i,j,k], :]
            tng = triangle.create_triang(r)
            # if scale R > reject_scale, reject
            if tng[5] < reject_scale:
                tl2.append(tng)

print len(tl1), len(tl2)
