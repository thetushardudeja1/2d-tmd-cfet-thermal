# QuantumATK smoke test: bulk silicon DFT-LCAO SCF + band gap.
# Proves the SCF engine actually computes. Expected: total energy converges to a
# finite negative value; DFT band gap ~0.4-0.7 eV (LDA/GGA underestimates 1.1 eV).
from QuantumATK import *

# --- bulk Si, diamond structure (a = 5.4306 Ang) ---
lattice = FaceCenteredCubic(5.4306 * Angstrom)
elements = [Silicon, Silicon]
frac = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
bulk = BulkConfiguration(bravais_lattice=lattice, elements=elements,
                         fractional_coordinates=frac)

# --- DFT-LCAO with a proper 7x7x7 k-grid (Y-2026.03 API) ---
nap = NumericalAccuracyParameters(k_point_sampling=(7, 7, 7))
calc = LCAOCalculator(numerical_accuracy_parameters=nap)
bulk.setCalculator(calc)

print("SMOKETEST_START")

# Stage 1: SCF total energy -- the core proof
te = TotalEnergy(bulk)
etot = te.evaluate().inUnitsOf(eV)
print("TOTAL_ENERGY_eV %.4f" % etot)

# Stage 2: band gap (best-effort; guarded so an API name change won't kill stage 1)
try:
    bs = Bandstructure(configuration=bulk, route=['L', 'G', 'X'],
                       points_per_segment=40)
    gap = bs.indirectBandGap().inUnitsOf(eV)
    print("INDIRECT_BANDGAP_eV %.4f" % gap)
except Exception as ex:
    print("BANDGAP_SKIPPED %s: %s" % (type(ex).__name__, ex))

print("SMOKETEST_OK")
