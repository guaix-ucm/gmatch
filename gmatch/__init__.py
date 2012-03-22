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

def gmatch(cat1, cat2, reject_scale=10.0, eps=1e-3):

    step = 0
    maxstep = 3
    ic1 = cat1
    ic2 = cat2
    while True:
        _logger.info('number of points cat1 %i', len(ic1))
        _logger.info('number of points cat2 %i', len(ic2))
        maxmatch = min(ic1.shape[0], ic2.shape[0])
        if maxmatch < 1:
            _logger.warning('no possible matches between catalogues')
            return None
        _logger.info('maximum number of matches %i', maxmatch)
        pm = gmatch_once(ic1, ic2, reject_scale, eps)

        nmatches = len(pm)
        _logger.info('matches: %i, maximum matches: %i', nmatches, maxmatch)
        if nmatches == maxmatch:
            _logger.info('objects in current input catalogues are all matched in %i steps', step)
            return ic1[pm[:,0]], ic2[pm[:,1]]

        elif nmatches == 0:
            _logger.info('no matches in catalogues')
            return None
        elif nmatches < maxmatch:
            if step < maxstep:
                _logger.info('number of matches < number of maximum matches')
                _logger.info('starting over with only matched points')             
                
                ic1 = cat1[pm[:,0]]
                ic2 = cat2[pm[:,1]]
                step += 1
            else:
                _logger.info('no convergence in matches after %i steps', maxstep)
                return None
                
        else:
            _logger.error('something is very wrong')
            return None


def gmatch_once(cat1s, cat2s, reject_scale=10.0, eps=1e-3):
    _logger.info('generating triangles in catalogue 1')
    tl1 = list(triangle.create_triang(cat1s, reject_scale=reject_scale, ep=eps))
    c = cat1s.shape[0]
    _logger.info('expected triangles %i', c * (c - 1) * (c - 2 ) / 6)
    _logger.info('created triangles %i', len(tl1))

    _logger.info('generating triangles in catalogue 2')
    tl2 = list(triangle.create_triang(cat2s, reject_scale=reject_scale, ep=eps))
    c = cat2s.shape[0]
    _logger.info('expected triangles %i', c * (c - 1) * (c - 2 ) / 6)
    _logger.info('created triangles %i', len(tl2))

    mrt1 = max(tl1, key=operator.attrgetter('tR'))
    mct1 = max(tl1, key=operator.attrgetter('tC'))
    _logger.debug('max R tolerance 1 %f', mrt1.tR)
    _logger.debug('max C tolerance 1 %f', mct1.tC)
    mrt2 = max(tl2, key=operator.attrgetter('tR'))
    mct2 = max(tl2, key=operator.attrgetter('tC'))
    _logger.debug('max R tolerance 2 %f', mrt2.tR)
    _logger.debug('max C tolerance 2 %f', mct2.tC)

    maxR = math.sqrt(mrt1.tR**2 + mrt2.tR**2)
    maxC = math.sqrt(mct1.tC**2 + mct2.tC**2)
    maxdis = math.sqrt(maxR**2 + maxC**2)
    _logger.info('max query tolerance in R space %f', maxR)
    _logger.info('max query tolerance in C space %f', maxC)
    _logger.info('max query tolerance in R-C space %f', maxdis)

    _logger.info('spliting R and C in catalogues')
    tspace1 = numpy.array([[tl.R, tl.C] for tl in tl1])    
    tspace2 = numpy.array([[tl.R, tl.C] for tl in tl2])

    _logger.info('finding closer triangles')
    # We use a kD tree to search in the R-C space
    _logger.debug('create kdtree...')
    kdtree = KDTree(tspace2)
    _logger.debug('done')

    _logger.debug('query in tree...')
    r = kdtree.query_ball_point(tspace1, r=maxdis)
    # r is an array of lists
    _logger.debug('done')
    
    matches1 = []
    _logger.info('checking matches')
    for i, l in numpy.ndenumerate(r):
        i0 = i[0]
        # triangle in first catalogue
        t1 = tl1[i0]
        # list of triangles in secodn catalogue
        t2s = [tl2[i1] for i1 in l]
        _logger.debug('trg %i in cat1 has %i neighbours', i0, len(t2s))
        # Return the closets trg in R-C space
        # or None
        mm = triangle.match_triangs(t1, t2s)
        if mm:
            matches1.append(mm)
    
    _logger.info('we have %i matches', len(matches1))
    if not matches1:
        _logger.info('no matches between the catalogues')
        return []
    
    # removing false matches
    _logger.info('filtering matches')
    matches = triangle.clean_matches(matches1)
    # voting to select matched coordinates
    _logger.info('voting')
    pm = triangle.votes(matches, cat1s.shape[0], cat2s.shape[0])    
    return pm

