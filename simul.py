
import numpy
import math

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

# Catalog1
idx = (50 < lista[:,0]) & (lista[:,0] < 150) & (50 < lista[:,1]) & (lista[:,1] < 150)

ll = lista[idx]
ll[:,0] -= 50
ll[:,1] -= 50

with open('cat1.txt', 'w') as fd:
    numpy.savetxt(fd, ll, fmt="%5.2f")

# Catalog2
idx = (0 < lista[:,0]) & (lista[:,0] < 100) & (0 < lista[:,1]) & (lista[:,1] < 100)
with open('cat2.txt', 'w') as fd:
    numpy.savetxt(fd, lista[idx], fmt="%5.2f")
