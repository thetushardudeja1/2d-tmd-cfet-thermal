"""PHASE 1d-B - is the melt-quenched HfO2 ACTUALLY amorphous (and physical)?

A too-short melt or too-slow quench can leave the structure crystalline; a bad quench can
leave it unphysically coordinated. Checks:
  1. DENSITY            target ~9.4 g/cm3 (amorphous HfO2 slightly below monoclinic 9.68)
  2. Hf-O BOND LENGTH   target ~2.1-2.2 A (monoclinic HfO2 Hf-O = 2.03-2.25 A)
  3. Hf-O COORDINATION  target ~7 (monoclinic is 7-fold; amorphous ~6.5-7.5)
  4. LONG-RANGE ORDER   compare the pair-distribution beyond ~6 A against the crystal:
                        crystal keeps sharp peaks, amorphous washes out  <-- the real test
"""
import numpy
from QuantumATK import *

CRYST = 'hfo2_crystal_ref.hdf5'


def pair_stats(conf, label):
    pos = numpy.array(conf.cartesianCoordinates().inUnitsOf(Angstrom))
    sym = [str(e.symbol()) for e in conf.elements()]
    cell = numpy.array(conf.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
    n = len(sym)

    # minimum-image distances (orthogonalised approximation is fine for these cells)
    Lx, Ly, Lz = cell[0][0], cell[1][1], cell[2][2]
    L = numpy.array([Lx, Ly, Lz])

    d_hfo, d_all = [], []
    for i in range(n):
        for j in range(i + 1, n):
            dv = pos[i] - pos[j]
            dv = dv - L * numpy.round(dv / L)
            r = float(numpy.sqrt((dv ** 2).sum()))
            if r > 12.0:
                continue
            d_all.append(r)
            if {sym[i], sym[j]} == {"Hf", "O"}:
                d_hfo.append(r)
    d_hfo = numpy.array(d_hfo)
    d_all = numpy.array(d_all)

    # 1. density
    vol = abs(numpy.linalg.det(cell)) * 1e-24                 # cm^3
    m_hf, m_o = 178.49, 16.00
    mass = (sym.count("Hf") * m_hf + sym.count("O") * m_o) / 6.02214076e23
    rho = mass / vol

    # 2/3. first Hf-O shell
    first = d_hfo[d_hfo < 2.8]
    n_hf = sym.count("Hf")
    coord = len(first) / max(n_hf, 1)

    print("--- %s ---" % label)
    print("  atoms %d (Hf %d, O %d)" % (n, sym.count("Hf"), sym.count("O")))
    print("  DENSITY        %.3f g/cm3      (target ~9.4)" % rho)
    print("  Hf-O mean      %.3f A          (target 2.1-2.2)" % (first.mean() if len(first) else -1))
    print("  Hf-O coord     %.2f            (target ~7)" % coord)

    # 4. long-range order: histogram beyond 6 A, report peak sharpness
    far = d_all[(d_all > 6.0) & (d_all < 12.0)]
    hist, edges = numpy.histogram(far, bins=60, range=(6.0, 12.0))
    hist = hist.astype(float)
    # normalise by shell volume so a random structure is flat
    centers = 0.5 * (edges[1:] + edges[:-1])
    hist /= centers ** 2
    hist /= hist.mean()
    sharpness = hist.std()          # crystal -> large; amorphous -> small
    print("  g(r>6A) STDDEV %.3f            (crystal >> amorphous)" % sharpness)
    return sharpness


print("VALIDATE_START")

# rebuild the crystalline reference for comparison
from SMW import *
mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (monoclinic)'](
    amorphous_density=9.390 * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2})
cryst = mat.pristineConfiguration().repeat(3, 3, 3)
s_cryst = pair_stats(cryst, "CRYSTALLINE monoclinic HfO2 (reference)")

amorph = nlread('amorphous_hfo2.hdf5', BulkConfiguration)[-1]
s_amorph = pair_stats(amorph, "MELT-QUENCHED HfO2")

print()
print("ORDER_RATIO  crystal/amorphous = %.2f" % (s_cryst / max(s_amorph, 1e-9)))
print("VERDICT: %s" % ("AMORPHOUS (long-range order destroyed)" if s_cryst / max(s_amorph, 1e-9) > 1.5
                       else "STILL ORDERED - melt longer / quench slower"))
print("VALIDATE_DONE")
