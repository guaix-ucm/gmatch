
from __future__ import print_function

import sys
import math
import itertools
import operator
import logging

import numpy
from scipy.spatial import cKDTree as KDTree

import triangle

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

def matching(cat1s, cat2s, reject_scale=10.0, eps=1e-3):
    common = min(cat1s.shape[0], cat2s.shape[0])
    logging.info('common number of points %i', common)

    print('generating triangles in catalogue 1')
    a = cat1s[:common,:]
    tl1 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

    print('generating triangles in catalogue 2')
    a = cat2s[:common,:]
    tl2 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

    print('expected triangles', common * (common - 1) * (common - 2 ) / 6)
    print('created triangles', 'cat1:', len(tl1), 'cat2:', len(tl2))

    mrt1 = max(tl1, key=operator.attrgetter('tR'))
    mct1 = max(tl1, key=operator.attrgetter('tC'))
    print('max R tolerance 1',mrt1.tR)
    print('max C tolerance 1',mct1.tC)
    mrt2 = max(tl2, key=operator.attrgetter('tR'))
    mct2 = max(tl2, key=operator.attrgetter('tC'))
    print('max R tolerance 2',mrt2.tR)
    print('max C tolerance 2',mct2.tC)

    maxR = math.sqrt(mrt1.tR**2 + mrt2.tR**2)
    maxC = math.sqrt(mct1.tC**2 + mct2.tC**2)
    maxdis = math.sqrt(maxR**2 + maxC**2)
    print('max query tolerance in R space', maxR)
    print('max query tolerance in C space', maxC)
    print('max query tolerance in R-C space', maxdis)

    print('spliting R and C in catalogues')
    tspace1 = numpy.array([[tl.R, tl.C] for tl in tl1])
    tspace2 = numpy.array([[tl.R, tl.C] for tl in tl2])

    print('create kdtree...', end=' ')
    kdtree = KDTree(tspace1)
    print('done')

    print('query in tree...', end=' ')
    r = kdtree.query(tspace2, distance_upper_bound=maxdis)
    print('done')

    _, indices = r

    print('checking matches')
    maxidx = len(tl1)
    matches1 = []
    for a,b in zip(tl2, indices):
        if b < maxidx:
            mt = triangle.match_triang(a, tl1[b])
            if mt:
                matches1.append(mt)

    print('filtering matches')
    matches = triangle.clean_matches(matches1)
    print('voting')
    pm = triangle.votes(matches, common)

    return pm

