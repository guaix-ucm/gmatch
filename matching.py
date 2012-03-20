
import sys
import math
import itertools
import operator
import logging

import numpy

from gmatch import matching

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

pm = matching(cat1s, cat2s)
    
for a,b in pm:
    print(cat1s[a], cat2s[b])

