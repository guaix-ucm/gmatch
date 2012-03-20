import sys
import math
import itertools
import operator

import numpy
from scipy.spatial import cKDTree as KDTree

import triangle

with open('cat1.txt') as fd:
    cat1 = numpy.loadtxt(fd)

with open('cat2.txt') as fd:
    cat2 = numpy.loadtxt(fd)

def sort_by_mag(a):
    # http://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
    return a[a[:,2].argsort()]


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

cat1s = sort_by_mag(cat1)
cat2s = sort_by_mag(cat2)
reject_scale = 10.0

common = min(cat1s.shape[0], cat2s.shape[0])
#common = 10
print 'Number of points', common

print 'generating triangles in catalogue 1'
a = cat1s[:common,:]
tl1 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

print 'generating triangles in catalogue 2'
a = cat2s[:common,:]
tl2 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

print 'expected', common * (common - 1) * (common - 2 ) / 6
print 'created', len(tl1), len(tl2)


mrt1 = max(tl1, key=operator.attrgetter('tR'))
mct1 = max(tl1, key=operator.attrgetter('tC'))
print 'max R tolerance 1',mrt1.tR
print 'max C tolerance 1',mct1.tC
mrt2 = max(tl2, key=operator.attrgetter('tR'))
mct2 = max(tl2, key=operator.attrgetter('tC'))
print 'max R tolerance 2',mrt2.tR
print 'max C tolerance 2',mct2.tC

maxdis = math.sqrt(mrt1.tR**2 + mrt2.tR**2 + mct1.tC**2 + mct2.tC**2)
print 'max query distance', maxdis

print 'spliting'
tspace1 = numpy.array([[tl.R, tl.C] for tl in tl1])
tspace2 = numpy.array([[tl.R, tl.C] for tl in tl2])

print 'create kdtree'
kdtree = KDTree(tspace1)
print 'done'

print 'query in tree'
r = kdtree.query(tspace2, distance_upper_bound=maxdis)
print 'done'

dis, indices = r

maxidx = len(tl1)
matches1 = []
print len(indices)
print indices
for a,b in zip(tl2, indices):
    if b < maxidx:
        mt = triangle.match_triang(a, tl1[b])
        if mt:
            matches1.append(mt)

matches = triangle.clean_matches(matches1)
