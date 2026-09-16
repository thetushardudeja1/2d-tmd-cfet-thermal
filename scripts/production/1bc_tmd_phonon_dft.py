# PHASE 1b / 1c - DFT-LCAO phonons for monolayer TMDs.
#   usage:  atkpython tmd_phonon_dft.py MoS2|WSe2
#
# 1b (MoS2): cross-check the SW force-field result (validated in Phase 1a) against DFT.
# 1c (WSe2): no shipped force field exists -> DFT is the defensible route.
#
# VALIDATION TARGETS (Gamma-point, monolayer, literature):
#   MoS2 : E"~286  E'~385  A1'~403  A2"~470   max ~477 cm^-1
#   WSe2 : E"~176  E'~249  A1'~250  A2"~308   max ~310 cm^-1
# The KEY physics for our project: WSe2's phonon spectrum is far SOFTER than MoS2's
# (lower Debye frequency, heavier atoms) -> that is the origin of its ultralow kappa.
import sys
from QuantumATK import *
import numpy

CM1 = 8065.54

mat = sys.argv[1] if len(sys.argv) > 1 else "MoS2"

# ---------------- structures: 1H monolayer ----------------
if mat == "MoS2":
    a, dz, els = 3.18, 1.5645, [Molybdenum, Sulfur, Sulfur]
elif mat == "WSe2":
    a, dz, els = 3.28, 1.67, [Tungsten, Selenium, Selenium]
else:
    raise SystemExit("unknown material %s" % mat)

c = 20.0
lattice = Hexagonal(a * Angstrom, c * Angstrom)
f = dz / c
frac = [[1.0/3.0, 2.0/3.0, 0.50],
        [2.0/3.0, 1.0/3.0, 0.50 + f],
        [2.0/3.0, 1.0/3.0, 0.50 - f]]
conf = BulkConfiguration(bravais_lattice=lattice, elements=els,
                         fractional_coordinates=frac)

print("DFT_PHONON_START %s" % mat)

# ---------------- DFT-LCAO calculator (2D: 1 k-point in vacuum direction) ----------------
nap = NumericalAccuracyParameters(k_point_sampling=(6, 6, 1))
calculator = LCAOCalculator(numerical_accuracy_parameters=nap)
conf.setCalculator(calculator)
conf.update()
print("SCF_OK  E_tot = %.4f eV" % TotalEnergy(conf).evaluate().inUnitsOf(eV))

# ---------------- relax (phonons need a force minimum) ----------------
try:
    conf = OptimizeGeometry(conf, max_forces=0.02*eV/Ang, max_steps=100,
                            optimize_cell=False)
    print("OPT_OK  E_tot = %.4f eV" % TotalEnergy(conf).evaluate().inUnitsOf(eV))
except Exception as ex:
    print("OPT_SKIPPED %s: %s" % (type(ex).__name__, ex))

# ---------------- dynamical matrix (3x3 in-plane supercell; enough for Gamma modes) ----
dm = DynamicalMatrix(
    configuration=conf,
    filename='%s_dft_phonon.hdf5' % mat,
    object_id='dm',
    calculator=calculator,
    repetitions=(3, 3, 1),
)
dm.update()
print("DM_OK")

# ---------------- band structure ----------------
pb = PhononBandstructure(dynamical_matrix=dm, route=['G', 'M', 'K', 'G'],
                         points_per_segment=30)
E = numpy.array(pb.evaluate().inUnitsOf(eV))
w = E * CM1
print("N_KPOINTS %d  N_BANDS %d" % (E.shape[0], E.shape[1]))
print("%s MAX_FREQ_cm-1  %.1f" % (mat, w.max()))
print("%s MIN_FREQ_cm-1  %.1f" % (mat, w.min()))
g = numpy.sort(w[0, :])
print("%s GAMMA_MODES_cm-1  %s" % (mat, " ".join("%.1f" % x for x in g)))

nlsave('%s_dft_phonon.hdf5' % mat, pb)
try:
    pdos = PhononDensityOfStates(dynamical_matrix=dm,
                                 qpoints=MonkhorstPackGrid(9, 9, 1))
    nlsave('%s_dft_phonon.hdf5' % mat, pdos)
    print("PDOS_OK")
except Exception as ex:
    print("PDOS_SKIPPED %s: %s" % (type(ex).__name__, ex))

print("DFT_PHONON_OK %s" % mat)
