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

def build(a1, a2, a3):
    print a1, a2, a3
    idx = ["12", "13", "23"]
    r = numpy.array([a1, a2, a3])

    dist = numpy.triu(distance.cdist(r, r))
    r = dist[dist>0]
    print r
    i1=r.argmax()
    print i1, idx[i1]
    r[i1] = 0
    i1=r.argmax()
    print i1, idx[i1]
    r[i1] = 0
    i1=r.argmax()
    print i1, idx[i1]
    r[i1] = 0
    sys.exit(0)

a = normaliza(cat1s)
l = cat1s.shape[0]
tl = []
for i in range(l):
    for j in range(i + 1, l):
        for k in range(j + 1, l):
            r = a[[i,j,k], :]
            tl.append(triangle.create_triang(r))

print len(tl)


#print normaliza(cat1s)
