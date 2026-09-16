"""Diagnose the melt-quenched glass: is it DESTROYED, or does it have a few bad pairs?

MINDIST alone is not enough. One pathological pair is a local defect that a short
minimisation fixes; hundreds means the structure is unusable. Measure the distribution.

Also checks whether MINDIST was a FALSE ALARM: my mindist() applies minimum-image
wrapping, and if MD left coordinates unwrapped by ~1e5 A, then dv - L*round(dv/L)
suffers catastrophic cancellation and reports garbage small distances.
"""
import sys
import numpy
from QuantumATK import *

fn = sys.argv[1] if len(sys.argv) > 1 else 'v8_glass.hdf5'
conf = nlread(fn, BulkConfiguration)[-1]
p = numpy.array(conf.cartesianCoordinates().inUnitsOf(Angstrom))
sym = [str(e.symbol()) for e in conf.elements()]
cl = numpy.array(conf.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
L = numpy.array([cl[0][0], cl[1][1], cl[2][2]])
n = len(sym)

print("FILE %s  atoms %d  cell %.2f x %.2f x %.2f" % (fn, n, L[0], L[1], L[2]))
print("RAW_COORD_RANGE x[%.1f,%.1f] y[%.1f,%.1f] z[%.1f,%.1f]"
      % (p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max(),
         p[:, 2].min(), p[:, 2].max()))
print("ANY_NAN %s   ANY_HUGE %s"
      % (bool(numpy.isnan(p).any()), bool((numpy.abs(p) > 1e4).any())))

# wrap cleanly first
f = (numpy.linalg.inv(cl.T) @ p.T).T
f -= numpy.floor(f)
p = (cl.T @ f.T).T

# full pair distance histogram (minimum image, now on WRAPPED coords)
bad, dists = [], []
for i in range(n):
    dv = p[i] - p
    for k in range(3):
        dv[:, k] -= L[k] * numpy.round(dv[:, k] / L[k])
    r = numpy.sqrt((dv ** 2).sum(axis=1))
    r[i] = 1e9
    dists.append(r[r < 6.0])
    for j in numpy.where(r < 1.5)[0]:
        if j > i:
            bad.append((i, j, sym[i], sym[j], float(r[j])))
d = numpy.concatenate(dists)

print()
print("PAIRS  <0.5A: %d   <1.0A: %d   <1.5A: %d   <2.0A: %d"
      % ((d < 0.5).sum() // 2, (d < 1.0).sum() // 2, (d < 1.5).sum() // 2, (d < 2.0).sum() // 2))
print("MIN_DISTANCE %.3f A" % d.min())
print()
print("WORST PAIRS (<1.5 A):")
for (i, j, a, b, r) in sorted(bad, key=lambda x: x[4])[:12]:
    print("   %4d(%2s) - %4d(%2s)  %.3f A" % (i, a, j, b, r))
print("   total bad pairs: %d out of %d atoms" % (len(bad), n))

# first-shell Hf-O to see if the network is otherwise sane
hf = [i for i, s in enumerate(sym) if s == "Hf"]
ox = [i for i, s in enumerate(sym) if s == "O"]
hfo = []
coord = []
for i in hf:
    dv = p[i] - p[ox]
    for k in range(3):
        dv[:, k] -= L[k] * numpy.round(dv[:, k] / L[k])
    r = numpy.sqrt((dv ** 2).sum(axis=1))
    hfo.append(r[r < 2.8])
    coord.append(int((r < 2.8).sum()))
hfo = numpy.concatenate(hfo)
coord = numpy.array(coord)
print()
print("Hf-O first shell  %.3f +- %.3f A   (target 2.1-2.2)" % (hfo.mean(), hfo.std()))
print("Hf coordination   %.2f +- %.2f     (target ~7)" % (coord.mean(), coord.std()))
print("VERDICT: %s" % ("A FEW LOCAL DEFECTS - fixable by minimisation"
                       if len(bad) < 0.02 * n else
                       "STRUCTURE DESTROYED - protocol change needed"))
print("DIAG_DONE")
