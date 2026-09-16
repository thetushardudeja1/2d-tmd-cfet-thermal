"""PHASE 1d-B (v8) - MoS2 / amorphous HfO2 in an ORTHORHOMBIC cell.

ROOT CAUSE OF v5-v7 (three failures, one cause):
  I kept taking fractional coordinates from a MONOCLINIC-derived cell and applying them
  to a HEXAGONAL cell. That is a SHEAR distortion, not a stretch. It produced 1.49 A
  contacts immediately (real HfO2: Hf-O 2.2, O-O 2.6) and the optimiser then drove them
  to 0.000 A. Fixing the hot anneal and the z-wrapping was correct but did not touch this.

FIX: use an ORTHORHOMBIC cell for everything. MoS2 has an exact orthorhombic
representation, so no shear is ever needed.
  hexagonal a1=(a,0), a2=(-a/2, a*sqrt3/2)  ->  ortho A1=a1=(a,0), A2=a1+2*a2=(0, a*sqrt3)
  ortho cell holds 2 formula units:
      Mo at (0, 1/3) and (1/2, 5/6)
      S  at (1/2, 1/6) and (0, 2/3)   (each doubled at +-dz)
  verified below: Mo-S must come out 2.41 A.

AMORPHOUS START: vendor `amorphousPrebuilder(configuration, density, sigma)` - it keeps
the input CELL SHAPE and sets the density, so we hand it an orthorhombic seed and never
shear anything.

usage: atkpython 1d_iface_v8_orthorhombic.py [nx] [ny] [melt_ps] [quench_ps] [gap]
"""
import sys
import numpy
from QuantumATK import *

NX = int(sys.argv[1]) if len(sys.argv) > 1 else 7      # MoS2 cells along a
NY = int(sys.argv[2]) if len(sys.argv) > 2 else 4      # MoS2 cells along b (= a*sqrt3)
MELT_PS = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
QUENCH_PS = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0
GAP = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0

A, DZ_S = 3.16, 1.5645
B = A * numpy.sqrt(3.0)
RHO, V_FU = 9.390, 210.49 / (6.02214076e23 * 9.390) * 1e24
VACUUM, DT = 18.0, 0.5
T_MELT = 3000.0   # HfO2 melts ~3100 K; 3500 K pushed the MTP out of its domain

LA, LB = NX * A, NY * B
AREA = LA * LB
print("V8_START  ortho cell %.3f x %.3f A  area %.1f A^2  dt %.2f fs" % (LA, LB, AREA, DT))

N_FU = int(round(AREA * 20.0 / V_FU))
LZ = N_FU * V_FU / AREA
print("HFO2  N_fu=%d (%d atoms)  Lz=%.3f A  rho=%.3f" % (N_FU, 3 * N_FU, LZ, RHO))

# ---------- 1. overlap-free orthorhombic seed on a grid ----------
n_at = 3 * N_FU
s = (AREA * LZ / n_at) ** (1.0 / 3.0)
nx, ny, nz = max(int(round(LA / s)), 1), max(int(round(LB / s)), 1), max(int(round(LZ / s)), 1)
while nx * ny * nz < n_at:
    nz += 1
sites = [(i / nx, j / ny, k / nz) for i in range(nx) for j in range(ny) for k in range(nz)]
rng = numpy.random.RandomState(1234)
rng.shuffle(sites)
sites = sites[:n_at]
els = [Hafnium] * N_FU + [Oxygen] * (2 * N_FU)
seed_cell = UnitCell([LA, 0, 0] * Angstrom, [0, LB, 0] * Angstrom, [0, 0, LZ] * Angstrom)
seed = BulkConfiguration(bravais_lattice=seed_cell, elements=els,
                         fractional_coordinates=numpy.array(sites))
print("SEED grid %dx%dx%d = %d sites, spacing ~%.2f A" % (nx, ny, nz, nx * ny * nz, s))


def mindist(conf, tag, periodic_z=True):
    p = numpy.array(conf.cartesianCoordinates().inUnitsOf(Angstrom))
    cl = numpy.array(conf.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
    Lv = numpy.array([cl[0][0], cl[1][1], cl[2][2]])
    d = 1e9
    for i in range(0, len(p), 2):
        dv = p[i] - p
        dv[:, 0] -= Lv[0] * numpy.round(dv[:, 0] / Lv[0])
        dv[:, 1] -= Lv[1] * numpy.round(dv[:, 1] / Lv[1])
        if periodic_z:
            dv[:, 2] -= Lv[2] * numpy.round(dv[:, 2] / Lv[2])
        r = numpy.sqrt((dv ** 2).sum(axis=1))
        r[r < 1e-6] = 1e9
        d = min(d, r.min())
    print("  MINDIST[%s] %.3f A" % (tag, d))
    return d


mindist(seed, "grid seed")

# ---------- 2. vendor amorphous prebuilder (keeps the ORTHORHOMBIC cell shape) ----------
try:
    pre = amorphousPrebuilder(seed, density=RHO * gram / cm ** 3, sigma=0.1)
    print("PREBUILD_OK atoms=%d" % len(pre.elements()))
except Exception as ex:
    print("PREBUILD_FAIL %s: %s -- falling back to jittered grid" % (type(ex).__name__, str(ex)[:160]))
    fr = numpy.array(sites) + rng.normal(0, 0.01, (n_at, 3))
    pre = BulkConfiguration(bravais_lattice=seed_cell, elements=els,
                            fractional_coordinates=fr - numpy.floor(fr))
mindist(pre, "prebuilt")

# v8 FAILED here: the melt-quench produced MINDIST 0.149 A. Machine-learned
# potentials are fitted near equilibrium and typically have NO SHORT-RANGE
# REPULSIVE WALL. At 3500 K the system explores geometries outside the training
# set, forces go flat/attractive, and atoms collapse into one another.
# FIX: merge a ZBL screened-nuclear-repulsion term to supply the physical 1/r core.
_mtp_set = QuantumATK_MTP_HfO2_2022()
try:
    _mtp_set.merge(ZBL_Tapered())
    print('ZBL_MERGED_OK')
except Exception as _ex:
    print('ZBL_MERGE_FAIL %s: %s' % (type(_ex).__name__, str(_ex)[:160]))
mtp = TremoloXCalculator(parameters=_mtp_set)
pre.setCalculator(mtp)
pre.update()

# ---------- 3. melt-quench at 0.5 fs ----------
iv = MaxwellBoltzmannDistribution(temperature=T_MELT * Kelvin,
                                  remove_center_of_mass_momentum=True,
                                  enforce_temperature=True)
melt = Langevin(time_step=DT * femtoSecond, reservoir_temperature=T_MELT * Kelvin,
                friction=0.05 * femtoSecond ** -1, initial_velocity=iv)
t1 = MolecularDynamics(pre, trajectory_filename='v8_melt.hdf5',
                       steps=int(MELT_PS * 1000 / DT), log_interval=4000, method=melt)
print("MELT_OK")
mindist(t1[-1], "after melt")   # localise: melt vs quench
quench = Langevin(time_step=DT * femtoSecond, reservoir_temperature=T_MELT * Kelvin,
                  friction=0.05 * femtoSecond ** -1,
                  heating_rate=-(T_MELT - 300.0) / QUENCH_PS * Kelvin / picoSecond)
t2 = MolecularDynamics(t1[-1], trajectory_filename='v8_quench.hdf5',
                       steps=int(QUENCH_PS * 1000 / DT), log_interval=4000, method=quench)
glass = t2[-1]
print("QUENCH_OK")

# wrap fully (still fully periodic here - no vacuum yet)
cl = numpy.array(glass.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
p = numpy.array(glass.cartesianCoordinates().inUnitsOf(Angstrom))
f = (numpy.linalg.inv(cl.T) @ p.T).T
f -= numpy.floor(f)
glass = BulkConfiguration(bravais_lattice=glass.bravaisLattice(),
                          elements=glass.elements(), fractional_coordinates=f)
mindist(glass, "glass")
nlsave('v8_glass.hdf5', glass)

# ---------- 4. add vacuum, keep slab where it is ----------
C = LZ + VACUUM
p = numpy.array(glass.cartesianCoordinates().inUnitsOf(Angstrom))
big = UnitCell([LA, 0, 0] * Angstrom, [0, LB, 0] * Angstrom, [0, 0, C] * Angstrom)
sub = BulkConfiguration(bravais_lattice=big, elements=list(glass.elements()),
                        cartesian_coordinates=p * Angstrom)
z_top = p[:, 2].max()
print("SUBSTRATE_Z min %.3f max %.3f  (cell C=%.2f, vacuum above)" % (p[:, 2].min(), z_top, C))

# ---------- 5. orthorhombic MoS2, with a bond-length self-check ----------
mo_f = [(0.0, 1.0 / 3.0), (0.5, 5.0 / 6.0)]
s_f = [(0.5, 1.0 / 6.0), (0.0, 2.0 / 3.0)]
els2 = list(sub.elements())
cart = [list(q) for q in p]
z_mo = z_top + GAP + DZ_S
mo_xy, s_xy = [], []
for i in range(NX):
    for j in range(NY):
        for (fx, fy) in mo_f:
            x, y = (i + fx) * A, (j + fy) * B
            mo_xy.append((x, y)); els2.append(Molybdenum); cart.append([x, y, z_mo])
        for (fx, fy) in s_f:
            x, y = (i + fx) * A, (j + fy) * B
            s_xy.append((x, y))
            els2 += [Sulfur, Sulfur]
            cart += [[x, y, z_mo + DZ_S], [x, y, z_mo - DZ_S]]

# self-check: nearest Mo-S must be 2.41 A
best = 1e9
for (mx, my) in mo_xy[:4]:
    for (sx, sy) in s_xy:
        dx, dy = mx - sx, my - sy
        dx -= LA * round(dx / LA); dy -= LB * round(dy / LB)
        best = min(best, (dx * dx + dy * dy + DZ_S ** 2) ** 0.5)
print("MOS2_BOND_CHECK  Mo-S = %.3f A  (must be 2.41)" % best)
if abs(best - 2.41) > 0.08:
    print("MOS2_LATTICE_WRONG - aborting")
    raise SystemExit

iface = BulkConfiguration(bravais_lattice=big, elements=els2,
                          cartesian_coordinates=numpy.array(cart) * Angstrom)
n_tot = len(els2)
print("IFACE_ATOMS %d (oxide %d + MoS2 %d)" % (n_tot, len(p), n_tot - len(p)))
if mindist(iface, "stacked", periodic_z=False) < 1.0:
    print("OVERLAP_ABORT")
    raise SystemExit

# ---------- 6. relax ----------
pot = TorchX_MACE_MPA_0_medium()
pot.merge(DispersionD3BJ())
iface.setCalculator(TremoloXCalculator(parameters=pot))
iface.update()
e0 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_BEFORE %.4f eV (%.4f eV/atom)" % (e0, e0 / n_tot))
iface = OptimizeGeometry(iface, max_forces=0.01 * eV / Ang, max_steps=500,
                         optimize_cell=False, trajectory_filename='v8_relax.hdf5')
e1 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_AFTER  %.4f eV (%.4f eV/atom)" % (e1, e1 / n_tot))

q = numpy.array(iface.cartesianCoordinates().inUnitsOf(Angstrom))
sy = [str(e.symbol()) for e in iface.elements()]
z_ox = max(q[i][2] for i, t in enumerate(sy) if t in ("Hf", "O"))
z_s = min(q[i][2] for i, t in enumerate(sy) if t == "S")
mz = [q[i][2] for i, t in enumerate(sy) if t in ("Mo", "S")]
print("FINAL_GAP   %.3f A  (target ~3)" % (z_s - z_ox))
print("MOS2_SPREAD %.3f A  (ideal %.3f)" % (max(mz) - min(mz), 2 * DZ_S))
nlsave('mos2_hfo2_v8.hdf5', iface)
print("V8_DONE")
