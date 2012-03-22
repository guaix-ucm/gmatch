
import sys
import math
import itertools
import operator
import logging

import numpy

sys.path.append("..") 

from gmatch import gmatch

logging.basicConfig(level=logging.INFO)

with open('cat1.txt') as fd:
    cat1 = numpy.loadtxt(fd)

with open('cat2.txt') as fd:
    cat2 = numpy.loadtxt(fd)

def sort_by_mag(a):
    # http://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
    return a[a[:,2].argsort()]


cat1s = sort_by_mag(cat1)
cat2s = sort_by_mag(cat2)
reject_scale = 10.0
eps = 1e-3

matches = gmatch(cat1s[:], cat2s[:], reject_scale=reject_scale, eps=eps)

if matches is not None:
    print matches[0]
    print matches[1]
