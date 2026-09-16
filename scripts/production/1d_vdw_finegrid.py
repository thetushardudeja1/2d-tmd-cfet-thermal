# PHASE 1d - DOES THE MLIP DESCRIBE THE vdW GAP CORRECTLY?
#
# Phonons being right does NOT guarantee interlayer binding is right, and TBC depends
# critically on interfacial coupling. Materials-Project-trained MLIPs are largely PBE,
# often WITHOUT dispersion -> the classic failure mode is: no binding, or an
# equilibrium gap that is far too large.
#
# TEST: 2H-MoS2 BILAYER interlayer binding curve E(d).
# REFERENCE (bulk 2H-MoS2):
#   * interlayer spacing (Mo-Mo, = c/2)  ~ 6.15 A   [a=3.16, c=12.29]
#   * interlayer binding energy          ~ 20-30 meV/atom
# A potential with NO vdW gives a very shallow/absent minimum at a much larger distance.
from QuantumATK import *
import numpy

A = 3.16          # in-plane lattice constant (Ang)
DZ = 1.5645       # Mo-S vertical offset within a layer
C = 40.0          # big vacuum so layers do not see periodic images
# FINE GRID (0.05 A) around the minimum, + far points for the asymptote
D_LIST = [5.6, 5.8, 6.0, 6.05, 6.10, 6.15, 6.20, 6.25, 6.30, 6.35, 6.40,
          6.45, 6.50, 6.55, 6.60, 6.65, 6.70, 6.80, 7.0, 7.5, 8.5, 12.0]

CANDIDATES = [
    ("MACE_MPA+D3",  "TorchX_MACE_MPA_0_medium",  True),
    ("MatterSim+D3", "TorchX_MatterSim_v1_0_0_5M", True),
]


def make_potential(pname, add_d3):
    pot = globals()[pname]()
    if add_d3:
        try:
            pot.merge(DispersionD3BJ())
        except Exception as ex:
            print("     (D3 merge failed: %s)" % ex)
    return pot


def bilayer(d):
    """2H-stacked MoS2 bilayer, Mo-Mo vertical separation = d."""
    lat = Hexagonal(A * Angstrom, C * Angstrom)
    z1 = 0.30 * C
    z2 = z1 + d
    els, cart = [], []
    # layer 1: Mo at (1/3,2/3)
    x1, y1 = A * (1/3. + 2/3. * 0.5), A * (2/3.) * (3 ** 0.5 / 2)
    # use fractional coords instead - simpler and unambiguous
    frac = []
    frac += [[1/3., 2/3., z1 / C], [2/3., 1/3., (z1 + DZ) / C], [2/3., 1/3., (z1 - DZ) / C]]
    # layer 2: 2H stacking -> Mo sits above S column
    frac += [[2/3., 1/3., z2 / C], [1/3., 2/3., (z2 + DZ) / C], [1/3., 2/3., (z2 - DZ) / C]]
    els = [Molybdenum, Sulfur, Sulfur, Molybdenum, Sulfur, Sulfur]
    return BulkConfiguration(bravais_lattice=lat, elements=els,
                             fractional_coordinates=frac)


def monolayer():
    lat = Hexagonal(A * Angstrom, C * Angstrom)
    z1 = 0.30 * C
    frac = [[1/3., 2/3., z1 / C], [2/3., 1/3., (z1 + DZ) / C], [2/3., 1/3., (z1 - DZ) / C]]
    return BulkConfiguration(bravais_lattice=lat,
                             elements=[Molybdenum, Sulfur, Sulfur],
                             fractional_coordinates=frac)


print("VDW_TEST_START")

for label, pname, add_d3 in CANDIDATES:
    print("=" * 74)
    print("POTENTIAL", label)
    try:
        # isolated monolayer reference
        m = monolayer()
        m.setCalculator(TremoloXCalculator(parameters=make_potential(pname, add_d3)))
        m.update()
        e_mono = TotalEnergy(m).evaluate().inUnitsOf(eV)

        curve = []
        for d in D_LIST:
            b = bilayer(d)
            b.setCalculator(TremoloXCalculator(parameters=make_potential(pname, add_d3)))
            b.update()
            e = TotalEnergy(b).evaluate().inUnitsOf(eV)
            # binding energy per atom, referenced to two isolated monolayers
            eb = (e - 2.0 * e_mono) / 6.0 * 1000.0     # meV/atom
            curve.append((d, eb))
            print("   d = %5.2f A   E_bind = %+8.2f meV/atom" % (d, eb))

        ds = numpy.array([c[0] for c in curve])
        es = numpy.array([c[1] for c in curve])
        i = int(numpy.argmin(es))
        print("   --> MINIMUM at d = %.2f A, E_bind = %.2f meV/atom" % (ds[i], es[i]))
        print("   --> REFERENCE   d ~ 6.15 A, E_bind ~ -20 to -30 meV/atom")
        dev = 100.0 * (ds[i] - 6.15) / 6.15
        verdict = "GOOD" if (abs(dev) < 8 and -45 < es[i] < -8) else "POOR (vdW not captured)"
        print("   --> d error %+.1f%%   VERDICT: %s" % (dev, verdict))
    except Exception as ex:
        print("   FAIL %s: %s" % (type(ex).__name__, str(ex)[:200]))

print("VDW_TEST_DONE")
