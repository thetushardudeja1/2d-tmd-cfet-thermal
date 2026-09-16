"""PHASE 1d-B - RIGOROUS amorphousness validation.

WHY v1 OF THIS TEST WAS BAD:
  It histogrammed pair distances over 6-12 A in a 3x3x3 cell that is only ~15.9 A wide.
  Half the box is ~8 A, so everything past ~8 A was dominated by PERIODIC IMAGES, which
  impose artificial order regardless of how amorphous the material is. The test could not
  pass at that cell size. Fixed here by (a) a bigger cell and (b) better metrics.

METRICS (all standard for amorphous-oxide characterisation):
  1. DENSITY                     -> 9.39 g/cm3 target
  2. Hf-O first-shell distance   -> 2.1-2.2 A
  3. COORDINATION DISTRIBUTION   -> crystal is a single sharp value (7); amorphous spreads
  4. O-Hf-O BOND-ANGLE SPREAD    -> the classic amorphous signature: sharp peaks broaden
  5. g(r) in a SAFE window (3-9 A, inside half-box) -> crystal sharp, amorphous washed out
  6. AMORPHISATION ENERGY vs the RELAXED crystal (not the strained one)

usage: atkpython 1d_validate_amorphous_v2.py <amorphous.hdf5> [nx ny nz]
"""
import sys
import numpy
from QuantumATK import *
from SMW import *

fname = sys.argv[1] if len(sys.argv) > 1 else 'amorphous_hfo2.hdf5'
NX = int(sys.argv[2]) if len(sys.argv) > 2 else 4
NY = int(sys.argv[3]) if len(sys.argv) > 3 else 4
NZ = int(sys.argv[4]) if len(sys.argv) > 4 else 4


def analyse(conf, label):
    pos = numpy.array(conf.cartesianCoordinates().inUnitsOf(Angstrom))
    sym = [str(e.symbol()) for e in conf.elements()]
    cell = numpy.array(conf.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
    L = numpy.array([cell[0][0], cell[1][1], cell[2][2]])
    n = len(sym)
    half = 0.5 * L.min()

    def mic(a, b):
        dv = pos[a] - pos[b]
        return dv - L * numpy.round(dv / L)

    hf = [i for i, s in enumerate(sym) if s == "Hf"]
    ox = [i for i, s in enumerate(sym) if s == "O"]

    # --- neighbour lists ---
    coord, angles, d_hfo, d_all = [], [], [], []
    for i in hf:
        nb = []
        for j in ox:
            dv = mic(i, j)
            r = float(numpy.sqrt((dv ** 2).sum()))
            if r < 2.8:
                nb.append(dv / r)
                d_hfo.append(r)
        coord.append(len(nb))
        for p in range(len(nb)):
            for q in range(p + 1, len(nb)):
                c = float(numpy.clip(numpy.dot(nb[p], nb[q]), -1, 1))
                angles.append(numpy.degrees(numpy.arccos(c)))
    for i in range(0, n, 2):                       # subsample for speed
        for j in range(i + 1, n, 2):
            dv = mic(i, j)
            r = float(numpy.sqrt((dv ** 2).sum()))
            if r < half:
                d_all.append(r)

    coord = numpy.array(coord, dtype=float)
    angles = numpy.array(angles)
    d_hfo = numpy.array(d_hfo)
    d_all = numpy.array(d_all)

    vol = abs(numpy.linalg.det(cell)) * 1e-24
    mass = (sym.count("Hf") * 178.49 + sym.count("O") * 16.00) / 6.02214076e23

    # g(r) in a SAFE window strictly inside half-box
    lo, hi = 3.0, min(9.0, half - 0.5)
    win = d_all[(d_all > lo) & (d_all < hi)]
    h, e = numpy.histogram(win, bins=50, range=(lo, hi))
    c = 0.5 * (e[1:] + e[:-1])
    h = h.astype(float) / c ** 2
    h /= h.mean()

    print("--- %s ---" % label)
    print("  atoms %d   box %.1f x %.1f x %.1f A   half-box %.1f A" % (n, L[0], L[1], L[2], half))
    print("  DENSITY           %.3f g/cm3" % (mass / vol))
    print("  Hf-O first shell  %.3f +- %.3f A" % (d_hfo.mean(), d_hfo.std()))
    print("  COORD mean/std    %.2f / %.2f   (crystal: sharp, std~0)" % (coord.mean(), coord.std()))
    print("  O-Hf-O angle std  %.2f deg      (crystal: sharp, amorphous: broad)" % angles.std())
    print("  g(r) %.1f-%.1f std  %.3f          (crystal >> amorphous)" % (lo, hi, h.std()))
    return {"coord_std": coord.std(), "angle_std": angles.std(), "gr_std": h.std(),
            "rho": mass / vol}


print("VALIDATE2_START")

mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (monoclinic)'](
    amorphous_density=9.390 * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2})
cryst = mat.pristineConfiguration().repeat(NX, NY, NZ)
c_stats = analyse(cryst, "CRYSTALLINE monoclinic HfO2 (relaxed density)")

amorph = nlread(fname, BulkConfiguration)[-1]
a_stats = analyse(amorph, "MELT-QUENCHED HfO2  [%s]" % fname)

# --- energies: amorphisation energy vs the RELAXED crystal ---
calc = TremoloXCalculator(parameters=QuantumATK_MTP_HfO2_2022())
cryst.setCalculator(calc); cryst.update()
e_c = TotalEnergy(cryst).evaluate().inUnitsOf(eV) / len(cryst.elements())
amorph.setCalculator(calc); amorph.update()
e_a = TotalEnergy(amorph).evaluate().inUnitsOf(eV) / len(amorph.elements())

print()
print("E_crystal(relaxed)  %.4f eV/atom" % e_c)
print("E_amorphous         %.4f eV/atom" % e_a)
print("AMORPHISATION_E     %+.4f eV/atom   (literature ~ +0.10 to +0.20)" % (e_a - e_c))

print()
print("SIGNATURES (amorphous should be LARGER for the first two, SMALLER for g(r)):")
print("  coord std   crystal %.3f -> amorphous %.3f" % (c_stats["coord_std"], a_stats["coord_std"]))
print("  angle std   crystal %.2f -> amorphous %.2f deg" % (c_stats["angle_std"], a_stats["angle_std"]))
print("  g(r) std    crystal %.3f -> amorphous %.3f" % (c_stats["gr_std"], a_stats["gr_std"]))

ok = (a_stats["angle_std"] > c_stats["angle_std"] * 1.15 and
      a_stats["gr_std"] < c_stats["gr_std"] * 0.7)
print("VERDICT: %s" % ("AMORPHOUS ✓" if ok else "NOT CONVINCINGLY AMORPHOUS - melt longer / quench slower / bigger cell"))
print("VALIDATE2_DONE")
