import math
import collections

import numpy as np

Triangle = collections.namedtuple('Triangle', ['v1', 'v2', 'v3', 'logp', 'hel', 'R', 'tR', 'C', 'tC'])

def norma(x):
    n = np.sqrt(np.dot(x, x.conj()))
    return n


def create_triang(v, ep=1e-3):
    # aristas
    a = v[[1,2,0], 0:2] - v[:, 0:2] # 1-0, 2-1, 0-2
    # normas de las aristas
    n = [norma(ar) for ar in a]
    # perimetro
    p = sum(n)
    ll = [(ni, ai, idx, (idx + 1) % 3) for ni, ai, idx in zip(n, a, range(3))]
    ls = sorted(ll)
    sides, aristas, idxs, nidxs = zip(*ls)

    ov = v[idxs, :]
    oa = np.array(aristas)

    e = np.cross(oa, [oa[1], oa[2], oa[0]])

    sg = np.sign(e)
    if np.any(sg != sg[0]):
        print 'reorder'
    R = sides[2] / sides[0]
    C = np.dot(oa[0], oa[2]) / (sides[2] * sides[0])
    dep1 = (1.0 / (sides[2])**2 + 1.0 / sides[0]**2 - C / (sides[2] * sides[0]))
    tR = 2 * R * R * ep * ep * dep1
    tC = 2 * (1 - C**2) * ep**2 * dep1 + 3 * C**2 * ep**4 * dep1**2

    return Triangle(v[0], v[1], v[2], math.log(p), sg[0], R, tR, C, tC)

if __name__ == '__main__':
    # vertice
    v = np.array([[0.16090669, 0.40127186,  15.25      ],
    [0.95456394,   0.90577636,  16.01      ],
    [0.21223826,   0.02564918,  16.14      ]])
    print create_triang(v, ep=1e-3)
