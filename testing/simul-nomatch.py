
import numpy
import math

numpy.random.seed(seed=232183)

coords = numpy.random.random(400)
coords.shape = (200, 2)
coords *= 200

class Expo(object):
    def __init__(self, m1, m2, al):
        self.alinv = 1.0 / al
        self.e1 = math.exp(al * m1)        
        self.e2 = math.exp(al * m2)        

    def __call__(self, x):
        return self.alinv * numpy.log(self.e1 + self.e2 * x)

e = Expo(12, 20, 1)

mags = numpy.random.random(200)
mags = e(mags)

lista = numpy.zeros((200, 3))
lista[:,0:2] = coords
lista[:,2] = mags

with open('master.txt', 'w') as fd:
    numpy.savetxt(fd, lista, fmt="%f")

# Catalog1
ll = []
for i in lista:
    if (0 < i[0] < 50) and (0 < i[1] < 50):
        ll.append(i)
ll = numpy.asarray(ll) - numpy.array([0, 0, 0])

with open('cat1.txt', 'w') as fd:
    numpy.savetxt(fd, ll, fmt="%f")

# Catalog2
coords = numpy.random.random(400)
coords.shape = (200, 2)
coords *= 200
lista[:,0:2] = coords
lista[:,2] = mags

ll = []
for i in lista:
    if (5 < i[0] < 55) and (5 < i[1] < 55):
        ll.append(i)
ll = numpy.asarray(ll) - numpy.array([5, 5, 0])
with open('cat2.txt', 'w') as fd:
    numpy.savetxt(fd, ll, fmt="%f")
