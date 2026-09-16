# PHASE 1a - VALIDATION GATE
# Monolayer 1H-MoS2 phonon band structure + DOS using the SHIPPED vendor potential
# StillingerWeber_MoS_2013, following the vendor template pattern from
#   _downloads/.../cnt_projected_phonon_dos.py   (DynamicalMatrix + TremoloXCalculator)
#
# VALIDATION TARGETS (published monolayer MoS2):
#   * max phonon frequency  ~473-490 cm^-1
#   * E' (E2g) ~385 cm^-1 , A1' (A1g) ~403 cm^-1   (Raman-active optical modes at Gamma)
#   * 3 acoustic branches -> ~0 at Gamma, ZA quadratic
#   * NO large imaginary (negative) modes => geometry is at equilibrium
from QuantumATK import *
import numpy

CM1 = 8065.54    # 1 eV in cm^-1

# ---------------- structure: monolayer 1H-MoS2 ----------------
a, c, dz = 3.18, 20.0, 1.5645          # Ang ; c = vacuum ; dz = Mo-S vertical offset
lattice = Hexagonal(a * Angstrom, c * Angstrom)
elements = [Molybdenum, Sulfur, Sulfur]
f = dz / c
frac = [[1.0/3.0, 2.0/3.0, 0.50],
        [2.0/3.0, 1.0/3.0, 0.50 + f],
        [2.0/3.0, 1.0/3.0, 0.50 - f]]
conf = BulkConfiguration(bravais_lattice=lattice, elements=elements,
                         fractional_coordinates=frac)

print("PHASE1A_START MoS2 StillingerWeber_MoS_2013")

# ---------------- calculator: shipped vendor force field ----------------
potentialSet = StillingerWeber_MoS_2013()
calculator = TremoloXCalculator(parameters=potentialSet)
conf.setCalculator(calculator)
conf.update()
print("CALC_OK  E_tot = %.4f eV" % TotalEnergy(conf).evaluate().inUnitsOf(eV))

# ---------------- relax to the potential's own equilibrium ----------------
# (phonons must be evaluated at a force minimum or acoustic modes go imaginary)
try:
    conf = OptimizeGeometry(conf, max_forces=0.01*eV/Ang, max_steps=200,
                            optimize_cell=False)
    print("OPT_OK   E_tot = %.4f eV" % TotalEnergy(conf).evaluate().inUnitsOf(eV))
except Exception as ex:
    print("OPT_SKIPPED %s: %s" % (type(ex).__name__, ex))

# ---------------- dynamical matrix (2D: repeat in-plane only) ----------------
dynamical_matrix = DynamicalMatrix(
    configuration=conf,
    filename='MoS2_phonon.hdf5',
    object_id='dm',
    calculator=calculator,
    repetitions=(5, 5, 1),
)
dynamical_matrix.update()
print("DM_OK")

# ---------------- phonon band structure ----------------
pb = PhononBandstructure(
    dynamical_matrix=dynamical_matrix,
    route=['G', 'M', 'K', 'G'],
    points_per_segment=40,
)
E = pb.evaluate().inUnitsOf(eV)          # (n_kpoints, n_bands)
E = numpy.array(E)
w = E * CM1

print("N_KPOINTS %d  N_BANDS %d" % (E.shape[0], E.shape[1]))
print("MAX_FREQ_cm-1     %.1f" % w.max())
print("MIN_FREQ_cm-1     %.1f   (large negative => imaginary modes)" % w.min())

# Gamma is the first point of the route
g = numpy.sort(w[0, :])
print("GAMMA_MODES_cm-1  " + " ".join("%.1f" % x for x in g))
print("GAMMA_OPTICAL_cm-1 " + " ".join("%.1f" % x for x in g[3:]))  # skip 3 acoustic

# ---------------- phonon DOS ----------------
try:
    pdos = PhononDensityOfStates(dynamical_matrix=dynamical_matrix,
                                 qpoints=MonkhorstPackGrid(11, 11, 1))
    print("PDOS_OK")
    nlsave('MoS2_phonon.hdf5', pdos)
except Exception as ex:
    print("PDOS_SKIPPED %s: %s" % (type(ex).__name__, ex))

nlsave('MoS2_phonon.hdf5', pb)
print("PHASE1A_OK")
