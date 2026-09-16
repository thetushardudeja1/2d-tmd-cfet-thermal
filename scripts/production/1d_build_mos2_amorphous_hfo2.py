"""PHASE 1d-B (final) - MoS2 on AMORPHOUS HfO2 interface.

KEY ADVANTAGE OF AMORPHOUS: it has no lattice, so we choose the cell to match MoS2
EXACTLY -> ZERO strain on the channel. (The crystalline route forced 0.54% strain and
a fragile Miller-index construction that failed.)

CELL CHOICE
  7 x MoS2 (a = 3.16 A) = 22.120 A in-plane, hexagonal.
  HfO2 slab is built into that SAME in-plane cell, with the out-of-plane length set so
  the density is exactly 9.390 g/cm3:
      V_per_formula_unit = 210.49 / (6.02214e23 * 9.390) cm3 = 37.217 A^3
      Lz = N_fu * 37.217 / A_inplane
  The 3500 K melt erases all memory of the (distorted) starting configuration, so
  reshaping the crystal into this cell before melting is legitimate.

PROTOCOL (D1): place MoS2 ~3 A above the slab, then relax to < 2e-2 eV/A.
POTENTIALS: MTP_HfO2 for the melt-quench (cheap); MACE_MPA + D3BJ for the interface
            relaxation (validated: phonons <1%, vdW binding -24.94 meV/atom).

usage: atkpython 1d_build_mos2_amorphous_hfo2.py [n_mos2] [melt_ps] [quench_ps] [gap_A]
"""
import sys
import numpy
from QuantumATK import *
from SMW import *

N_MOS2 = int(sys.argv[1]) if len(sys.argv) > 1 else 7
MELT_PS = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
QUENCH_PS = float(sys.argv[3]) if len(sys.argv) > 3 else 43.0
GAP = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0

A_MOS2, DZ_S = 3.16, 1.5645
RHO = 9.390
V_FU = 210.49 / (6.02214076e23 * RHO) * 1e24      # A^3 per HfO2 formula unit
VACUUM = 15.0

L = N_MOS2 * A_MOS2
# hexagonal in-plane cell area
A_IN = L * L * numpy.sqrt(3.0) / 2.0
print("IFACE_START  N_MoS2=%d  L=%.3f A  area=%.1f A^2" % (N_MOS2, L, A_IN))

# ---------- 1. HfO2 formula units for a ~14 A slab at the right density ----------
LZ_TARGET = 14.0
N_FU = int(round(A_IN * LZ_TARGET / V_FU))
LZ = N_FU * V_FU / A_IN
print("HFO2  N_fu=%d (%d atoms)  Lz=%.3f A  -> rho=%.3f g/cm3"
      % (N_FU, 3 * N_FU, LZ, N_FU * 210.49 / (6.02214076e23 * A_IN * LZ * 1e-24)))

# ---------- 2. seed from the correct crystal, reshaped into our cell ----------
mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (monoclinic)'](
    amorphous_density=RHO * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2})
seed = mat.pristineConfiguration()
nrep = int(numpy.ceil((N_FU / 4.0) ** (1.0 / 3.0)))      # 4 f.u. per monoclinic cell
seed = seed.repeat(nrep, nrep, nrep)
sy = [str(e.symbol()) for e in seed.elements()]
fr = numpy.array(seed.fractionalCoordinates())
hf_i = [i for i, s in enumerate(sy) if s == "Hf"]
o_i = [i for i, s in enumerate(sy) if s == "O"]
# take exactly N_FU Hf and 2*N_FU O
keep = hf_i[:N_FU] + o_i[:2 * N_FU]
els = [Hafnium] * N_FU + [Oxygen] * (2 * N_FU)
frac = fr[keep]
frac = frac - numpy.floor(frac)

hexcell = UnitCell([L, 0.0, 0.0] * Angstrom,
                   [-L / 2.0, L * numpy.sqrt(3.0) / 2.0, 0.0] * Angstrom,
                   [0.0, 0.0, LZ] * Angstrom)
slab = BulkConfiguration(bravais_lattice=hexcell, elements=els,
                         fractional_coordinates=frac)
print("SEED_ATOMS %d" % len(els))

mtp = TremoloXCalculator(parameters=QuantumATK_MTP_HfO2_2022())
slab.setCalculator(mtp); slab.update()

# ---------- 3. melt-quench (erases the reshaping distortion) ----------
iv = MaxwellBoltzmannDistribution(temperature=3500 * Kelvin,
                                  remove_center_of_mass_momentum=True,
                                  enforce_temperature=True)
melt = Langevin(time_step=1 * femtoSecond, reservoir_temperature=3500 * Kelvin,
                friction=0.05 * femtoSecond ** -1, initial_velocity=iv)
t1 = MolecularDynamics(slab, trajectory_filename='if_melt.hdf5',
                       steps=int(MELT_PS * 1000), log_interval=2000, method=melt)
molten = t1[-1]
print("MELT_OK")

rate = -(3500.0 - 300.0) / QUENCH_PS
quench = Langevin(time_step=1 * femtoSecond, reservoir_temperature=3500 * Kelvin,
                  friction=0.05 * femtoSecond ** -1,
                  heating_rate=rate * Kelvin / picoSecond)
t2 = MolecularDynamics(molten, trajectory_filename='if_quench.hdf5',
                       steps=int(QUENCH_PS * 1000), log_interval=2000, method=quench)
amorph = t2[-1]
print("QUENCH_OK")
nlsave('amorphous_hfo2_iface.hdf5', amorph)

# ---------- 4. stack MoS2 on top, with vacuum ----------
pos = numpy.array(amorph.cartesianCoordinates().inUnitsOf(Angstrom))
print("RAW_SLAB_Z  min %.3f  max %.3f  (cell Lz %.3f)" % (pos[:,2].min(), pos[:,2].max(), LZ))

# CRITICAL: MD trajectories store UNWRAPPED coordinates. After 83 ps the slab has drifted
# far outside [0, Lz]; using z.max() directly put MoS2 INSIDE the oxide (v1 of this script:
# E_before +138 eV/atom, gap 0.904 A, MoS2 spread 11.1 A). Wrap in ALL directions first.
_cellv = numpy.array(amorph.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
_inv = numpy.linalg.inv(_cellv.T)
_f = (_inv @ pos.T).T
_f = _f - numpy.floor(_f)                 # wrap into [0,1)
pos = (_cellv.T @ _f.T).T
print("WRAPPED_SLAB_Z  min %.3f  max %.3f" % (pos[:,2].min(), pos[:,2].max()))
z_top = pos[:, 2].max()
C = LZ + GAP + 2 * DZ_S + VACUUM
els2 = list(amorph.elements())
cart = [list(p) for p in pos]

w1 = numpy.array([A_MOS2, 0.0])
w2 = numpy.array([-A_MOS2 / 2.0, A_MOS2 * numpy.sqrt(3.0) / 2.0])
z_mo = z_top + GAP + DZ_S
for i in range(N_MOS2):
    for j in range(N_MOS2):
        pm = (i + 1.0 / 3.0) * w1 + (j + 2.0 / 3.0) * w2
        ps = (i + 2.0 / 3.0) * w1 + (j + 1.0 / 3.0) * w2
        els2 += [Molybdenum, Sulfur, Sulfur]
        cart += [[pm[0], pm[1], z_mo],
                 [ps[0], ps[1], z_mo + DZ_S],
                 [ps[0], ps[1], z_mo - DZ_S]]

cell2 = UnitCell([L, 0.0, 0.0] * Angstrom,
                 [-L / 2.0, L * numpy.sqrt(3.0) / 2.0, 0.0] * Angstrom,
                 [0.0, 0.0, C] * Angstrom)
iface = BulkConfiguration(bravais_lattice=cell2, elements=els2,
                          cartesian_coordinates=numpy.array(cart) * Angstrom)
n_tot = len(els2)
print("IFACE_ATOMS %d  (HfO2 %d + MoS2 %d)  cell height %.2f A"
      % (n_tot, 3 * N_FU, 3 * N_MOS2 * N_MOS2, C))

# ---------- 5. relax with the production potential ----------
pot = TorchX_MACE_MPA_0_medium()
pot.merge(DispersionD3BJ())
mace = TremoloXCalculator(parameters=pot)
# --- FAIL FAST: nearest-neighbour check. Any pair < 1.0 A means the stack is broken
# and relaxing it just produces garbage (as in v1 of this script).
_p = numpy.array(iface.cartesianCoordinates().inUnitsOf(Angstrom))
_L = numpy.array([L, L * numpy.sqrt(3.0) / 2.0, C])
_dmin, _npairs = 1e9, 0
for _i in range(0, len(_p), 3):
    _dv = _p[_i] - _p
    _dv[:, 2] -= C * numpy.round(_dv[:, 2] / C)
    _r = numpy.sqrt((_dv ** 2).sum(axis=1))
    _r[_r < 1e-6] = 1e9
    _dmin = min(_dmin, _r.min())
    _npairs += int((_r < 1.2).sum())
print("MIN_DISTANCE %.3f A   pairs<1.2A: %d" % (_dmin, _npairs))
if _dmin < 1.0:
    print("OVERLAP_ABORT - structure is broken, not relaxing")
    raise SystemExit

iface.setCalculator(mace); iface.update()
e0 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_BEFORE %.4f eV (%.4f eV/atom)" % (e0, e0 / n_tot))

iface = OptimizeGeometry(iface, max_forces=0.02 * eV / Ang, max_steps=400,
                         optimize_cell=False,
                         trajectory_filename='iface_relax.hdf5')
e1 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_AFTER  %.4f eV (%.4f eV/atom)   drop %.4f eV" % (e1, e1 / n_tot, e1 - e0))

# ---------- 6. geometry report ----------
p = numpy.array(iface.cartesianCoordinates().inUnitsOf(Angstrom))
s = [str(e.symbol()) for e in iface.elements()]
z_ox = max(p[i][2] for i, q in enumerate(s) if q in ("Hf", "O"))
z_s = min(p[i][2] for i, q in enumerate(s) if q == "S")
mz = [p[i][2] for i, q in enumerate(s) if q in ("Mo", "S")]
print("FINAL_GAP     %.3f A  (oxide top -> lowest S)" % (z_s - z_ox))
print("MOS2_SPREAD   %.3f A  (ideal %.3f)" % (max(mz) - min(mz), 2 * DZ_S))

nlsave('mos2_amorphous_hfo2.hdf5', iface)
print("SAVED mos2_amorphous_hfo2.hdf5")
print("IFACE_DONE")
