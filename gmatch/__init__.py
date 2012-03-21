
import sys
import math
import itertools
import operator
import logging

import numpy
# cKDTree does not supoort
# query ball
from scipy.spatial import KDTree

import triangle

def normalize1(a):
    mi = a.min()
    mx = a.max()

    a = (a - mi) / (mx - mi)
    return a

def normalize(a):
    b = a.copy()
    b[:,0] = normalize1(b[:,0])
    b[:,1] = normalize1(b[:,1])

    return b

def matching(cat1s, cat2s, nmatch=None, reject_scale=10.0, eps=1e-3):
    common = min(cat1s.shape[0], cat2s.shape[0])
    if nmatch is not None:
        common = min(nmatch, common)

    logging.info('common number of points %i', common)

    logging.info('generating triangles in catalogue 1')
    a = normalize(cat1s[:common,:])
    tl1 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

    logging.info('generating triangles in catalogue 2')
    a = normalize(cat2s[:common,:])
    tl2 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

    logging.info('expected triangles %i', common * (common - 1) * (common - 2 ) / 6)
    logging.info('created triangles', 'cat1: %i cat2: %i', len(tl1), len(tl2))

    mrt1 = max(tl1, key=operator.attrgetter('tR'))
    mct1 = max(tl1, key=operator.attrgetter('tC'))
    logging.info('max R tolerance 1 %f', mrt1.tR)
    logging.info('max C tolerance 1 %f', mct1.tC)
    mrt2 = max(tl2, key=operator.attrgetter('tR'))
    mct2 = max(tl2, key=operator.attrgetter('tC'))
    logging.info('max R tolerance 2 %f', mrt2.tR)
    logging.info('max C tolerance 2 %f',mct2.tC)

    maxR = math.sqrt(mrt1.tR**2 + mrt2.tR**2)
    maxC = math.sqrt(mct1.tC**2 + mct2.tC**2)
    maxdis = math.sqrt(maxR**2 + maxC**2)
    logging.info('max query tolerance in R space', maxR)
    logging.info('max query tolerance in C space', maxC)
    logging.info('max query tolerance in R-C space', maxdis)

    logging.info('spliting R and C in catalogues')
    tspace1 = numpy.array([[tl.R, tl.C] for tl in tl1])
    tspace2 = numpy.array([[tl.R, tl.C] for tl in tl2])

    for c,v in zip(tl1, tl2):
        logging.info(c,v)

    logging.info('create kdtree...', end=' ')
    kdtree = KDTree(tspace1)
    logging.info('done')

    logging.info('query in tree...', end=' ')
    r = kdtree.query_ball_point(tspace2, r=maxdis)
    # r is an array of lists
    logging.info('done')
    
    matches1 = []
    logging.info('checking matches')
    for i, l in numpy.ndenumerate(r):
        i0 = i[0]
        t2 = tl2[i0]
        t1s = [tl1[i1] for i1 in l]
        mm = triangle.match_triangs(t2, t1s)
        if mm:
            matches1.append(mm)
    
    logging.info('we have', len(matches1),'matches')
    logging.info('filtering matches')
    matches = triangle.clean_matches(matches1)
    logging.info('voting')
    pm = triangle.votes(matches, common)

    if len(pm) < common:
        logging.info('we should start over')

    return pm

