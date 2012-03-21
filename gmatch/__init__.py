#
# Copyright 2012 Universidad Complutense de Madrid
# 
# This file is part of Gmatch
# 
# Gmatch is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Gmatch is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Gmatch.  If not, see <http://www.gnu.org/licenses/>.
#

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

_logger = logging.getLogger('gmatch')

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

    _logger.info('common number of points %i', common)

    _logger.info('generating triangles in catalogue 1')
    a = normalize(cat1s[:common,:])
    tl1 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

    _logger.info('generating triangles in catalogue 2')
    a = normalize(cat2s[:common,:])
    tl2 = list(triangle.create_triang(a, common, reject_scale=reject_scale))

    _logger.info('expected triangles %i', common * (common - 1) * (common - 2 ) / 6)
    _logger.info('created triangles', 'cat1: %i cat2: %i', len(tl1), len(tl2))

    mrt1 = max(tl1, key=operator.attrgetter('tR'))
    mct1 = max(tl1, key=operator.attrgetter('tC'))
    _logger.info('max R tolerance 1 %f', mrt1.tR)
    _logger.info('max C tolerance 1 %f', mct1.tC)
    mrt2 = max(tl2, key=operator.attrgetter('tR'))
    mct2 = max(tl2, key=operator.attrgetter('tC'))
    _logger.info('max R tolerance 2 %f', mrt2.tR)
    _logger.info('max C tolerance 2 %f',mct2.tC)

    maxR = math.sqrt(mrt1.tR**2 + mrt2.tR**2)
    maxC = math.sqrt(mct1.tC**2 + mct2.tC**2)
    maxdis = math.sqrt(maxR**2 + maxC**2)
    _logger.info('max query tolerance in R space', maxR)
    _logger.info('max query tolerance in C space', maxC)
    _logger.info('max query tolerance in R-C space', maxdis)

    _logger.info('spliting R and C in catalogues')
    tspace1 = numpy.array([[tl.R, tl.C] for tl in tl1])
    tspace2 = numpy.array([[tl.R, tl.C] for tl in tl2])

    if _logger.isEnabledFor(logging.DEBUG):
        for c,v in zip(tl1, tl2):
            _logger.debug('%s %s',c,v)

    _logger.info('create kdtree...')
    kdtree = KDTree(tspace1)
    _logger.info('done')

    _logger.info('query in tree...')
    r = kdtree.query_ball_point(tspace2, r=maxdis)
    # r is an array of lists
    _logger.info('done')
    
    matches1 = []
    _logger.info('checking matches')
    for i, l in numpy.ndenumerate(r):
        i0 = i[0]
        t2 = tl2[i0]
        t1s = [tl1[i1] for i1 in l]
        mm = triangle.match_triangs(t2, t1s)
        if mm:
            matches1.append(mm)
    
    _logger.info('we have %i matches', len(matches1))
    _logger.info('filtering matches')
    matches = triangle.clean_matches(matches1)
    _logger.info('voting')
    pm = triangle.votes(matches, common)

    if len(pm) < common:
        _logger.info('we should start over')

    return pm

