import sys
import math
import itertools

import numpy
import scipy.spatial.distance as distance

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

a = cat1s[:common,:]
l = common
tl1 = []
for idx in itertools.combinations(range(common), 3):
    r = a[idx, :]
    tng = triangle.create_triang(r)
    # if scale R > reject_scale, reject
    if tng[5] < reject_scale:
        tl1.append(tng)
    else:
        print 'reject'

for tl in tl1:
    print tl

#a = normaliza(cat2s[:common,:])
a = cat2s[:common,:]
l = common
tl2 = []
for idx in itertools.combinations(range(common), 3):
    r = a[idx, :]
    tng = triangle.create_triang(r)
    # if scale R > reject_scale, reject
    if tng[5] < reject_scale:
        tl2.append(tng)
    else:
        print 'reject'

print 'expected', common * (common - 1) * (common - 2 ) / 6
print len(tl1), len(tl2)

first = tl1[0]

def naive_match(tg1, tg2):
    R1 = tg1[5]
    tR1 = tg1[6]
    C1 = tg1[7]
    tC1 = tg1[8]
    R2 = tg2[5]
    tR2 = tg2[6]
    C2 = tg2[7]
    tC2 = tg2[8]

    sen1 = (R1 - R2)**2 - tR1**2 - tR2**2
    # print 'sen1', sen1

    if sen1 > 0:
        return None

    sen2 = (C1 - C2)**2 - tC1**2 - tC2**2

    if sen2 > 0:
        return None
    # print 'sen2', sen2

    return sen1, sen2


for first in tl1:
    for tl in tl2:
        m = naive_match(first, tl)
        if m is not None:
            print 't1=',first
            print 't2=',tl
