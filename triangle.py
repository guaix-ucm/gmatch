import math
import collections
import itertools

import numpy as np

Triangle = collections.namedtuple('Triangle', ['v0', 'v1', 'v2', 'i0', 'i1', 'i2', 'logp', 'hel', 'R', 'tR', 'C', 'tC'])

MatchedTriangles = collections.namedtuple('MatchedTriangles', ['t0', 't1', 'hel', 'logm'])

def norma(x):
    n = np.sqrt(np.dot(x, x.conj()))
    return n

def votes(matches, catsize):

    # shape of the catalogues, not of the mactches
    nm = catsize

    vot = np.zeros((nm, nm), dtype='int')
    pair = np.empty((nm,2), dtype='int')
    pair.fill(-1)

    result = []

    for m in matches:
        t0 = m.t0
        t1 = m.t1

        vot[t0.i0, t1.i0] += 1
        vot[t0.i1, t1.i1] += 1
        vot[t0.i2, t1.i2] += 1

    vmx = vot.max()
    if vmx <= 0:
        print 'no match'
    
    sortv = np.argsort(vot, axis=None)
    id0, id1 = np.unravel_index(sortv[::-1], (nm, nm))
    for i,j in zip(id0, id1):
        val = vot[i,j]
        if val <= 0:
            # votes are 0
            print 'votes are 0, ending'
            break
            
        if 2 * val < vmx:
            print 'votes are a half of the maximum, ending'
            # votes are a half of the maximum
            break
        if pair[i,0] != -1 or pair[j,1] != -1:
            # the point is already matched
            print 'point',i,j,'already matched, ending'
            break

        pair[i,0] = j
        pair[j,1] = i
        result.append((i, j))
    return result

def _scale_factor(mf, mt):
    scale = 0
    if mf > mt:
        scale = 1
    elif 0.1 * mt > mf:
        scale = 3
    else:
        scale = 2
    return scale

def clean_matches(matches):

    nmatches = len(matches)
    nnmatches = nmatches

    while True:
        npl = nm = 0
        logm = []
        for match in matches:
            if match.hel > 0:
                npl += 1
            elif match.hel < 0:
                nm -= 1
            else:
                break
            logm.append(match.logm)

        print 'n+ is', npl
        print 'n- is', nm

        mt = abs(npl - nm)
        mf = npl + nm - mt

        scale = _scale_factor(mf, mt)
        print 'scale factor is', scale

        lgmrr = np.array(logm)

        med = lgmrr.mean()
        std = lgmrr.std()

        print 'log M, average=',med, ' std=',std

        if std == 0:
            print 'std is 0, end matching'
            break

        newmatches = []
        print 'removing false matches due to scale'
        for match in matches:
            z = (match.logm - med ) / (scale * std)
            if -1 <= z <= 1:
                newmatches.append(match)

        print 'matches were', nmatches
        nnmatches = len(newmatches)
        print 'matches are', nnmatches

        if nmatches == nnmatches:
            matches = newmatches
            break
        matches = newmatches


    return matches

def match_triang(t1, t2):

    def mr(t1, t2):
        return (t1.R - t2.R)**2 - t1.tR**2 - t2.tR**2
    def mc(t1, t2):
        return (t1.C - t2.C)**2 - t1.tC**2 - t2.tC**2

    sen1 = mr(t1, t2)

    if sen1 > 0:
        return None

    sen2 = mc(t1, t2)

    if sen2 > 0:
        return None

    return MatchedTriangles(t1, t2, t1.hel * t2.hel, t1.logp - t2.logp)
    

def create_triang(vlist, many, reject_scale=10, ep=1e-3):
    for idx in itertools.combinations(range(many), 3):
        t = create_triang_(vlist, idx, ep)
        if t.R < reject_scale:
            yield t

def create_triang_(vlist, idx, ep=1e-3):
    v = vlist[idx, :]
    # aristas
    a = v[[1,2,0], 0:2] - v[:, 0:2] # 1-0, 2-1, 0-2
    # normas de las aristas
    n = [norma(ar) for ar in a]
    # perimetro
    p = sum(n)
    ll = [(ni, ai, ids, (ids + 1) % 3) for ni, ai, ids in zip(n, a, range(3))]
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

    return Triangle(v[0], v[1], v[2], idx[0], idx[1], idx[2], math.log(p), sg[0], R, tR, C, tC)

if __name__ == '__main__':
    # vertice
    v = np.array([[0.16090669, 0.40127186,  15.25      ],
    [0.95456394,   0.90577636,  16.01      ],
    [0.21223826,   0.02564918,  16.14      ]])
    print create_triang(v, ep=1e-3)
