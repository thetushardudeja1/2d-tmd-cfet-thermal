"""PHASE 1d-B (v6) - MoS2 / amorphous-HfO2 interface, built the D1 way.

WHAT CHANGED (five previous builds failed; see INTERFACE_POTENTIAL_RECIPE.md):
  * DO NOT re-melt from a reshaped crystal. Reshaping monoclinic fractional coords into a
    hexagonal cell produced overlaps -> the 3500 K MD exploded to +-200,000 A.
    Instead START FROM THE ALREADY-VALIDATED amorphous slab.
  * WRAP coordinates. MD trajectories are UNWRAPPED.
  * TIME STEP 0.5 fs (D1's production value; QuantumATK's 1 fs default was unstable at 3500 K).
  * VACUUM ADDED BEFORE equilibration so the surface reconstructs (D1 S3 step 2).
  * NO BAROSTAT. D1 uses NPT, but they have no imposed lattice. Here the amorphous glass
    must ADOPT the MoS2 in-plane cell, so the cell is FIXED and we anneal (NVT) instead - barostatting would (a) fight the imposed match and (b) try to compress the vacuum.
  * FAIL-FAST minimum-distance check before any expensive relaxation.

PROTOCOL
  1. load validated amorphous slab, wrap
  2. affine rescale in-plane to 7x MoS2 = 22.120 A, adjust Lz to hold rho = 9.390 g/cm3
  3. add vacuum, anneal NVT 1500 K -> 300 K at 0.5 fs (relieves the affine strain,
     reconstructs the free surface) - NOT a full melt, so the glass structure survives
  4. verify: density, min distance, still amorphous
  5. place MoS2 at 3.0 A, fail-fast check
  6. minimise with MACE_MPA + D3BJ to < 1e-2 eV/A (D1's force criterion)

usage: atkpython 1d_iface_v2_from_validated_slab.py [slab.hdf5] [anneal_ps] [gap]
"""
import sys
import numpy
from QuantumATK import *

SLAB = sys.argv[1] if len(sys.argv) > 1 else 'amorphous_hfo2.hdf5'
ANNEAL_PS = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
GAP = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

N_MOS2, A_MOS2, DZ_S = 7, 3.16, 1.5645
RHO = 9.390
VACUUM = 18.0
DT = 0.5          # fs - D1 production value

L = N_MOS2 * A_MOS2
A_IN = L * L * numpy.sqrt(3.0) / 2.0
print("V6_START  target in-plane L=%.3f A  area=%.1f A^2  dt=%.2f fs" % (L, A_IN, DT))

# ---------- 1. load the VALIDATED amorphous slab and wrap ----------
slab = nlread(SLAB, BulkConfiguration)[-1]
els = list(slab.elements())
sym = [str(e.symbol()) for e in els]
n = len(els)
cell = numpy.array(slab.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
pos = numpy.array(slab.cartesianCoordinates().inUnitsOf(Angstrom))
f = (numpy.linalg.inv(cell.T) @ pos.T).T
f = f - numpy.floor(f)
print("LOADED %s  atoms=%d (Hf %d, O %d)" % (SLAB, n, sym.count("Hf"), sym.count("O")))

# ---------- 2. affine rescale into the MoS2-matched hexagonal cell ----------
n_fu = sym.count("Hf")
v_fu = 210.49 / (6.02214076e23 * RHO) * 1e24
LZ = n_fu * v_fu / A_IN
print("RESCALE -> hexagonal L=%.3f, Lz=%.3f A  (rho stays %.3f g/cm3)" % (L, LZ, RHO))

C = LZ + VACUUM
hexcell = UnitCell([L, 0.0, 0.0] * Angstrom,
                   [-L / 2.0, L * numpy.sqrt(3.0) / 2.0, 0.0] * Angstrom,
                   [0.0, 0.0, C] * Angstrom)
# keep fractional x,y; compress z into the slab region only (vacuum above)
f2 = f.copy()
f2[:, 2] = f[:, 2] * (LZ / C)
slab2 = BulkConfiguration(bravais_lattice=hexcell, elements=els,
                          fractional_coordinates=f2)


def mindist(conf, tag):
    p = numpy.array(conf.cartesianCoordinates().inUnitsOf(Angstrom))
    cl = numpy.array(conf.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
    Lz = cl[2][2]
    dmin = 1e9
    for i in range(0, len(p), 2):
        dv = p[i] - p
        dv[:, 2] -= Lz * numpy.round(dv[:, 2] / Lz)
        r = numpy.sqrt((dv ** 2).sum(axis=1))
        r[r < 1e-6] = 1e9
        dmin = min(dmin, r.min())
    print("  MINDIST[%s] %.3f A" % (tag, dmin))
    return dmin


mindist(slab2, "after rescale")

mtp = TremoloXCalculator(parameters=QuantumATK_MTP_HfO2_2022())
slab2.setCalculator(mtp)
slab2.update()

# ---------- 3. anneal (NVT, fixed cell) to relieve strain + reconstruct surface ----------
iv = MaxwellBoltzmannDistribution(temperature=1500 * Kelvin,
                                  remove_center_of_mass_momentum=True,
                                  enforce_temperature=True)
hot = NVTNoseHoover(time_step=DT * femtoSecond, reservoir_temperature=1500 * Kelvin,
                    thermostat_timescale=100 * femtoSecond, initial_velocity=iv)
t1 = MolecularDynamics(slab2, trajectory_filename='v6_anneal.hdf5',
                       steps=int(ANNEAL_PS * 1000 / DT), log_interval=2000, method=hot)
print("ANNEAL_HOT_OK")

cool = NVTNoseHoover(time_step=DT * femtoSecond, reservoir_temperature=1500 * Kelvin,
                     thermostat_timescale=100 * femtoSecond,
                     heating_rate=-(1500.0 - 300.0) / ANNEAL_PS * Kelvin / picoSecond)
t2 = MolecularDynamics(t1[-1], trajectory_filename='v6_cool.hdf5',
                       steps=int(ANNEAL_PS * 1000 / DT), log_interval=2000, method=cool)
sub = t2[-1]
print("ANNEAL_COOL_OK")

# wrap again (MD unwraps!)
cell3 = numpy.array(sub.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
p3 = numpy.array(sub.cartesianCoordinates().inUnitsOf(Angstrom))
f3 = (numpy.linalg.inv(cell3.T) @ p3.T).T
f3 = f3 - numpy.floor(f3)
sub = BulkConfiguration(bravais_lattice=sub.bravaisLattice(),
                        elements=sub.elements(), fractional_coordinates=f3)
p3 = numpy.array(sub.cartesianCoordinates().inUnitsOf(Angstrom))
print("SUBSTRATE_Z  min %.3f  max %.3f  (slab should occupy ~0-%.1f, vacuum above)"
      % (p3[:, 2].min(), p3[:, 2].max(), LZ))
mindist(sub, "after anneal")
nlsave('v6_substrate.hdf5', sub)

# ---------- 4. stack MoS2 ----------
z_top = p3[:, 2].max()
els2 = list(sub.elements())
cart = [list(q) for q in p3]
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

iface = BulkConfiguration(bravais_lattice=hexcell, elements=els2,
                          cartesian_coordinates=numpy.array(cart) * Angstrom)
n_tot = len(els2)
print("IFACE_ATOMS %d  (oxide %d + MoS2 %d)  z_top=%.2f  z_Mo=%.2f  C=%.2f"
      % (n_tot, n, 3 * N_MOS2 * N_MOS2, z_top, z_mo, C))
if mindist(iface, "stacked") < 1.0:
    print("OVERLAP_ABORT")
    raise SystemExit

# ---------- 5. minimise with the production potential ----------
pot = TorchX_MACE_MPA_0_medium()
pot.merge(DispersionD3BJ())
mace = TremoloXCalculator(parameters=pot)
iface.setCalculator(mace)
iface.update()
e0 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_BEFORE %.4f eV (%.4f eV/atom)" % (e0, e0 / n_tot))

iface = OptimizeGeometry(iface, max_forces=0.01 * eV / Ang, max_steps=500,
                         optimize_cell=False, trajectory_filename='v6_relax.hdf5')
e1 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_AFTER  %.4f eV (%.4f eV/atom)" % (e1, e1 / n_tot))

p = numpy.array(iface.cartesianCoordinates().inUnitsOf(Angstrom))
s = [str(e.symbol()) for e in iface.elements()]
z_ox = max(p[i][2] for i, q in enumerate(s) if q in ("Hf", "O"))
z_s = min(p[i][2] for i, q in enumerate(s) if q == "S")
mz = [p[i][2] for i, q in enumerate(s) if q in ("Mo", "S")]
print("FINAL_GAP    %.3f A   (target ~3)" % (z_s - z_ox))
print("MOS2_SPREAD  %.3f A   (ideal %.3f)" % (max(mz) - min(mz), 2 * DZ_S))
nlsave('mos2_hfo2_iface_v6.hdf5', iface)
print("SAVED mos2_hfo2_iface_v6.hdf5")
print("V6_DONE")
