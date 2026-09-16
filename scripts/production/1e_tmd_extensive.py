# EXTENSIVE TMD PHONON STUDY
#   usage: atkpython tmd_extensive.py <MATERIAL> [REPS]
#
# Systematic first-principles phonon study across the TMD family relevant to CFET
# channel selection (n-candidates: MoS2, WS2, MoSe2 ; p-candidates: WSe2, MoTe2, WTe2).
#
# Extracts, per material:
#   * Gamma optical modes (E", E', A1', A2")      -> literature validation
#   * max phonon frequency (Debye-like scale)      -> phonon energy scale
#   * LA / TA acoustic SOUND VELOCITIES near Gamma -> the quantity that sets kappa
#     (kappa ~ (1/3) C v^2 tau, so v is the first-principles driver of the kappa ranking)
#   * phonon DOS
#
# Sanity targets: MoS2 v_LA ~6.6 km/s, v_TA ~4.2 km/s ; WSe2 v_LA ~3.3 km/s
import sys
from QuantumATK import *
import numpy

CM1 = 8065.54          # 1 eV -> cm^-1
EV2RADS = 1.51927e15   # 1 eV -> rad/s  (E/hbar)

# 1H monolayer geometries: a (in-plane, Ang), dz (metal-chalcogen vertical offset, Ang)
MATS = {
    "MoS2":  (3.18, 1.5645, [Molybdenum, Sulfur,   Sulfur]),
    "WS2":   (3.18, 1.57,   [Tungsten,   Sulfur,   Sulfur]),
    "MoSe2": (3.32, 1.67,   [Molybdenum, Selenium, Selenium]),
    "WSe2":  (3.28, 1.67,   [Tungsten,   Selenium, Selenium]),
    "MoTe2": (3.55, 1.82,   [Molybdenum, Tellurium, Tellurium]),
    "WTe2":  (3.55, 1.80,   [Tungsten,   Tellurium, Tellurium]),
}

mat = sys.argv[1]
reps = int(sys.argv[2]) if len(sys.argv) > 2 else 3
a, dz, els = MATS[mat]
c = 20.0

lattice = Hexagonal(a * Angstrom, c * Angstrom)
f = dz / c
conf = BulkConfiguration(
    bravais_lattice=lattice, elements=els,
    fractional_coordinates=[[1/3., 2/3., 0.50],
                            [2/3., 1/3., 0.50 + f],
                            [2/3., 1/3., 0.50 - f]])

print("EXT_START %s reps=%d a=%.3f" % (mat, reps, a))

nap = NumericalAccuracyParameters(k_point_sampling=(6, 6, 1))
calculator = LCAOCalculator(numerical_accuracy_parameters=nap)
conf.setCalculator(calculator)
conf.update()
print("%s SCF_OK E=%.4f eV" % (mat, TotalEnergy(conf).evaluate().inUnitsOf(eV)))

try:
    conf = OptimizeGeometry(conf, max_forces=0.02*eV/Ang, max_steps=100, optimize_cell=False)
    print("%s OPT_OK E=%.4f eV" % (mat, TotalEnergy(conf).evaluate().inUnitsOf(eV)))
except Exception as ex:
    print("%s OPT_SKIPPED %s" % (mat, ex))

dm = DynamicalMatrix(configuration=conf, filename='%s_ext_r%d.hdf5' % (mat, reps),
                     object_id='dm', calculator=calculator, repetitions=(reps, reps, 1))
dm.update()
print("%s DM_OK" % mat)

NPTS = 60
pb = PhononBandstructure(dynamical_matrix=dm, route=['G', 'M', 'K', 'G'],
                         points_per_segment=NPTS)
E = numpy.array(pb.evaluate().inUnitsOf(eV))     # (nk, nbands)
w = E * CM1

print("%s MAXFREQ_cm-1 %.1f" % (mat, w.max()))
print("%s MINFREQ_cm-1 %.1f" % (mat, w.min()))
g = numpy.sort(w[0, :])
print("%s GAMMA_cm-1 %s" % (mat, " ".join("%.1f" % x for x in g)))

# ---- acoustic sound velocities from the Gamma->M segment slope ----
# hexagonal: |Gamma->M| = 2*pi/(a*sqrt(3))
a_m = a * 1e-10
dk = (2.0 * numpy.pi / (a_m * numpy.sqrt(3.0))) / NPTS      # 1/m per index step
Es = numpy.sort(E, axis=1)                                   # sort bands at each k
fitn = 6                                                     # points near Gamma for the fit
idx = numpy.arange(1, fitn + 1)
for b, label in [(0, "ZA"), (1, "TA"), (2, "LA")]:
    omega = Es[1:fitn+1, b] * EV2RADS                        # rad/s
    kk = idx * dk
    slope = numpy.polyfit(kk, omega, 1)[0]                   # m/s
    print("%s V_%s_m_s %.1f  (%.2f km/s)" % (mat, label, slope, slope/1000.0))

nlsave('%s_ext_r%d.hdf5' % (mat, reps), pb)
try:
    pdos = PhononDensityOfStates(dynamical_matrix=dm, qpoints=MonkhorstPackGrid(9, 9, 1))
    nlsave('%s_ext_r%d.hdf5' % (mat, reps), pdos)
    print("%s PDOS_OK" % mat)
except Exception as ex:
    print("%s PDOS_SKIPPED %s" % (mat, ex))

print("EXT_OK %s" % mat)
