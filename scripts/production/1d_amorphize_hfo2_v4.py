"""PHASE 1d-B (v3) - amorphous HfO2 by explicit melt-quench MD.

WHY THIS ROUTE
  v2 failed: AmorphousLayer is a *pure specification object* (only amorphousDensity /
  length / materialSpecifications / method / stoichiometry) - it builds nothing itself.
  It must be handed to MultilayerBuilder. But MultilayerBuilder wants a
  per-INTERFACE MTP (the vendor uses MTP_Si_SiO2, MTP_HfO2_SiO2, MTP_HfO2_TiN) and no
  MoS2/HfO2 MTP exists. So we do the melt-quench ourselves, transparently.

  Crystal structure comes from the SMW database (correct crystallography, NOT hand-built - v1 failed exactly by hand-building fluorite).

MELT-QUENCH follows the vendor HKMG example's Langevin annealing pattern:
  stage 1: Langevin @ melt T          (melt)
  stage 2: Langevin with heating_rate<0 down to 300 K   (quench)

CALCULATOR: QuantumATK_MTP_HfO2_2022 (purpose-trained for HfO2, cheap).

usage: atkpython 1d_amorphize_hfo2_md.py [nx] [ny] [nz] [melt_ps] [quench_ps]
"""
import sys
from QuantumATK import *
from SMW import *

NX = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NY = int(sys.argv[2]) if len(sys.argv) > 2 else 3
NZ = int(sys.argv[3]) if len(sys.argv) > 3 else 3
MELT_PS = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
QUENCH_PS = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
T_MELT = 3500.0

print("AMD_START %dx%dx%d melt=%.1fps quench=%.1fps" % (NX, NY, NZ, MELT_PS, QUENCH_PS))

# ---------- 1. correct monoclinic HfO2 from the SMW database ----------
mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (monoclinic)'](
    amorphous_density=9.390 * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2},
)
bulk = None
for acc in ["pristineConfiguration", "pristine_configuration", "configuration",
            "bulkConfiguration"]:
    if hasattr(mat, acc):
        try:
            out = getattr(mat, acc)
            out = out() if callable(out) else out
            if out is not None:
                bulk = out
                print("PRISTINE_VIA %s" % acc)
                break
        except Exception as ex:
            print("  %s failed: %s" % (acc, str(ex)[:120]))
if bulk is None:
    print("NO_PRISTINE - methods:", [m for m in dir(mat) if not m.startswith("_")][:30])
    raise SystemExit

n0 = len(bulk.elements())
print("PRISTINE_ATOMS %d" % n0)
try:
    L = bulk.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom)
    print("PRISTINE_CELL  a=%.3f b=%.3f c=%.3f A"
          % (sum(x*x for x in L[0])**0.5, sum(x*x for x in L[1])**0.5,
             sum(x*x for x in L[2])**0.5))
except Exception:
    pass

# ---------- 2. supercell ----------
sc = bulk.repeat(NX, NY, NZ)
n = len(sc.elements())
print("SUPERCELL_ATOMS %d" % n)

# ---------- 2b. EXPAND to the amorphous density (CRITICAL) ----------
# v3 FAILED because it melted at the CRYSTAL density (10.109 g/cm3). Amorphous HfO2 is
# LESS dense (9.390 g/cm3, the vendor HKMG value) - a compressed cell cannot disorder,
# and the validation showed long-range order surviving (ORDER_RATIO only 1.17).
import numpy as _np
_cell = _np.array(sc.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
_vol = abs(_np.linalg.det(_cell)) * 1e-24
_sym = [str(e.symbol()) for e in sc.elements()]
_mass = (_sym.count("Hf") * 178.49 + _sym.count("O") * 16.00) / 6.02214076e23
_rho0 = _mass / _vol
_scale = (_rho0 / 9.390) ** (1.0 / 3.0)
print("DENSITY_CRYSTAL %.3f g/cm3 -> target 9.390 ; linear scale %.5f" % (_rho0, _scale))
_frac = sc.fractionalCoordinates()
sc = BulkConfiguration(
    bravais_lattice=UnitCell(_cell[0]*_scale*Angstrom,
                             _cell[1]*_scale*Angstrom,
                             _cell[2]*_scale*Angstrom),
    elements=sc.elements(),
    fractional_coordinates=_frac)
_cell2 = _np.array(sc.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
print("DENSITY_EXPANDED %.3f g/cm3" % (_mass / (abs(_np.linalg.det(_cell2)) * 1e-24)))

calc = TremoloXCalculator(parameters=QuantumATK_MTP_HfO2_2022())
sc.setCalculator(calc)
sc.update()
print("E_CRYSTAL %.4f eV (%.4f eV/atom)" % (
    TotalEnergy(sc).evaluate().inUnitsOf(eV),
    TotalEnergy(sc).evaluate().inUnitsOf(eV) / n))

# ---------- 3. melt ----------
iv = MaxwellBoltzmannDistribution(temperature=T_MELT * Kelvin,
                                  remove_center_of_mass_momentum=True,
                                  enforce_temperature=True)
melt = Langevin(time_step=1 * femtoSecond,
                reservoir_temperature=T_MELT * Kelvin,
                friction=0.05 * femtoSecond ** -1,
                initial_velocity=iv)
print("MELT_RUN %d steps @ %.0f K" % (int(MELT_PS * 1000), T_MELT))
traj = MolecularDynamics(sc, trajectory_filename='hfo2_melt.hdf5',
                         steps=int(MELT_PS * 1000), log_interval=1000, method=melt)
molten = traj[-1] if hasattr(traj, "__getitem__") else traj
print("MELT_OK")

# ---------- 4. quench ----------
rate = -(T_MELT - 300.0) / QUENCH_PS
quench = Langevin(time_step=1 * femtoSecond,
                  reservoir_temperature=T_MELT * Kelvin,
                  friction=0.05 * femtoSecond ** -1,
                  heating_rate=rate * Kelvin / picoSecond)
print("QUENCH_RUN %d steps, rate %.1f K/ps" % (int(QUENCH_PS * 1000), rate))
traj2 = MolecularDynamics(molten, trajectory_filename='hfo2_quench.hdf5',
                          steps=int(QUENCH_PS * 1000), log_interval=1000, method=quench)
amorph = traj2[-1] if hasattr(traj2, "__getitem__") else traj2
print("QUENCH_OK")

# ---------- 5. relax + report ----------
amorph.setCalculator(calc)
amorph.update()
e = TotalEnergy(amorph).evaluate().inUnitsOf(eV)
print("E_AMORPHOUS %.4f eV (%.4f eV/atom)" % (e, e / n))
nlsave('amorphous_hfo2.hdf5', amorph)
print("SAVED amorphous_hfo2.hdf5")
print("AMD_DONE")
