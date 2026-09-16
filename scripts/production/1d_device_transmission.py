"""PHASE 1d-C - device phonon transmission -> TBC for TMD/HfO2 interfaces.

WHY THE STRUCTURE IS REBUILT HERE
  The relaxed slabs (FINAL_iface_*) are monolayer-on-substrate WITH VACUUM. Landauer
  transport needs SEMI-INFINITE LEADS on both sides, and phonons cannot propagate through
  vacuum. So the device is built fully periodic along z:
        [ bulk HfO2 electrode | HfO2 + vdW gap + TMD | bulk TMD electrode ]
  This is legitimate and directly comparable to D1: their Fig. 3 shows TBC RISES with TMD
  layer count and SATURATES beyond 3 layers at G = 13.0 +- 0.8 MW/m^2/K. A bulk-TMD lead is
  therefore the saturated multilayer limit, i.e. D1's saturated value - not an approximation.

KEY API (TOOLCHAIN_VALIDATION.md 4a/4b):
  * PhononTransmissionSpectrum takes `dynamical_matrix` FIRST (vendor example is stale)
  * DynamicalMatrix(use_equivalent_bulk=False) is REQUIRED for a device

usage: atkpython 1d_device_transmission.py MoS2|WSe2 [n_ox] [n_tmd] [gap]
"""
import sys
import numpy
from QuantumATK import *

TMD = sys.argv[1] if len(sys.argv) > 1 else 'MoS2'
N_OX = int(sys.argv[2]) if len(sys.argv) > 2 else 3      # HfO2 [111] repeats in central region
N_TMD = int(sys.argv[3]) if len(sys.argv) > 3 else 3     # TMD layers in central region
GAP = float(sys.argv[4]) if len(sys.argv) > 4 else None  # default: use relaxed value

_P = {'MoS2': (3.16, 1.5645, Molybdenum, Sulfur,   6.15, 4.697),
      'WSe2': (3.28, 1.67,   Tungsten,   Selenium, 6.49, 4.193)}
A_TMD, DZ, M_AT, X_AT, C_TMD, GAP_REL = _P[TMD]
if GAP is None:
    GAP = GAP_REL

print("DEVTX_START %s  n_ox=%d n_tmd=%d gap=%.3f" % (TMD, N_OX, N_TMD, GAP))

# ---------- 1. HfO2(111), O-terminated, NO vacuum (electrode-compatible) ----------
from SMW import *
mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (cubic)'](
    amorphous_density=9.390 * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2})
ox_unit = mat.pristineConfiguration().cleave(1, 1, 1, layers=1, outer_atom_index=1)
ov = numpy.array(ox_unit.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
a_ox, c_ox = numpy.linalg.norm(ov[0]), ov[2][2]
print("OX_UNIT atoms=%d a=%.4f  c_repeat=%.4f" % (len(ox_unit.elements()), a_ox, c_ox))

# ---------- 2. in-plane lattice match ----------
best = None
for m in range(1, 14):
    for n in range(1, 16):
        mis = abs(n * A_TMD - m * a_ox) / (m * a_ox)
        if mis < 0.02 and (best is None or mis < best[0]):
            best = (mis, n, m)
mis, N_T, M_O = best
L = M_O * a_ox
print("MATCH %d x %s (%.3f) vs %d x HfO2 (%.3f) -> %.3f%%  L=%.3f"
      % (N_T, TMD, N_T * A_TMD, M_O, L, 100 * mis, L))

v1 = ov[0] * M_O
v2 = ov[1] * M_O

# ---------- 3. assemble: oxide block (bottom) + gap + TMD block (top), periodic in z ----
els, cart = [], []
oxp = numpy.array(ox_unit.cartesianCoordinates().inUnitsOf(Angstrom))
oxs = list(ox_unit.elements())
n_ox_rep = N_OX + 2                       # +2 -> the electrode repeats at the bottom
for k in range(n_ox_rep):
    for i in range(M_O):
        for j in range(M_O):
            sh = i * ov[0] + j * ov[1] + numpy.array([0, 0, k * c_ox])
            for e, q in zip(oxs, oxp):
                els.append(e); cart.append(list(q + sh))
z_ox_top = max(c[2] for c in cart)

# TMD: 2H stacking, N_TMD+2 layers (+2 -> electrode repeats at the top)
w1, w2 = v1 / N_T, v2 / N_T
n_tmd_rep = N_TMD + 2
d_tmd = C_TMD / 2.0                       # interlayer spacing
z0 = z_ox_top + GAP + DZ
for k in range(n_tmd_rep):
    zc = z0 + k * d_tmd
    off = (w1 + w2) / 3.0 if k % 2 else numpy.zeros(3)   # 2H: alternate registry
    for i in range(N_T):
        for j in range(N_T):
            pm = (i + 1.0/3.0) * w1 + (j + 2.0/3.0) * w2 + off
            ps = (i + 2.0/3.0) * w1 + (j + 1.0/3.0) * w2 + off
            els += [M_AT, X_AT, X_AT]
            cart += [[pm[0], pm[1], zc], [ps[0], ps[1], zc + DZ], [ps[0], ps[1], zc - DZ]]

C = max(c[2] for c in cart) + DZ + d_tmd - DZ            # close the periodic stack
cell = UnitCell(v1 * Angstrom, v2 * Angstrom, [0.0, 0.0, C] * Angstrom)
central = BulkConfiguration(bravais_lattice=cell, elements=els,
                            cartesian_coordinates=numpy.array(cart) * Angstrom)
n_tot = len(els)
print("CENTRAL atoms=%d  C=%.3f  z_ox_top=%.3f" % (n_tot, C, z_ox_top))

# sanity
p = numpy.array(central.cartesianCoordinates().inUnitsOf(Angstrom))
dmin = 1e9
for i in range(0, len(p), 3):
    dv = p[i] - p
    dv[:, 0] -= v1[0] * numpy.round(dv[:, 0] / v1[0])
    dv[:, 2] -= C * numpy.round(dv[:, 2] / C)
    r = numpy.sqrt((dv ** 2).sum(axis=1)); r[r < 1e-6] = 1e9
    dmin = min(dmin, r.min())
print("MIN_DISTANCE %.3f A" % dmin)
if dmin < 1.0:
    print("SANITY_ABORT"); raise SystemExit

pot = TorchX_MACE_MPA_0_medium(); pot.merge(DispersionD3BJ())
calc = TremoloXCalculator(parameters=pot)
central.setCalculator(calc); central.update()
print("E_central %.4f eV (%.4f eV/atom)"
      % (TotalEnergy(central).evaluate().inUnitsOf(eV),
         TotalEnergy(central).evaluate().inUnitsOf(eV) / n_tot))

# ---------- 4. electrodes ----------
ox_elec = ox_unit.repeat(M_O, M_O, 2)
tmd_cell = UnitCell(v1 * Angstrom, v2 * Angstrom, [0.0, 0.0, C_TMD] * Angstrom)
te, tc = [], []
for k in range(2):
    zc = k * d_tmd
    off = (w1 + w2) / 3.0 if k % 2 else numpy.zeros(3)
    for i in range(N_T):
        for j in range(N_T):
            pm = (i + 1.0/3.0) * w1 + (j + 2.0/3.0) * w2 + off
            ps = (i + 2.0/3.0) * w1 + (j + 1.0/3.0) * w2 + off
            te += [M_AT, X_AT, X_AT]
            tc += [[pm[0], pm[1], zc], [ps[0], ps[1], zc + DZ], [ps[0], ps[1], zc - DZ]]
tmd_elec = BulkConfiguration(bravais_lattice=tmd_cell, elements=te,
                             cartesian_coordinates=numpy.array(tc) * Angstrom)
print("ELECTRODES  oxide=%d atoms   %s=%d atoms"
      % (len(ox_elec.elements()), TMD, len(tmd_elec.elements())))

dev = DeviceConfiguration(central, [ox_elec, tmd_elec])
dev.setCalculator(calc); dev.update()
print("DEVICE_OK")

# ---------- 5. dynamical matrix + transmission ----------
dm = DynamicalMatrix(configuration=dev, filename='tx_%s.hdf5' % TMD, object_id='dm',
                     calculator=calc, repetitions=(1, 1, 1),
                     use_equivalent_bulk=False,      # REQUIRED for a device
                     constrain_electrodes=False)
dm.update()
print("DM_OK")

pts = PhononTransmissionSpectrum(
    dynamical_matrix=dm,
    energies=numpy.linspace(0, 0.08, 161) * eV,      # 0-645 cm-1 covers both TMD + HfO2
    qpoints=MonkhorstPackGrid(5, 5),
    infinitesimal=1e-6 * eV,
    self_energy_calculator=RecursionSelfEnergy())
nlsave('tx_%s.hdf5' % TMD, pts)
print("PTS_OK")

# ---------- 6. Landauer TBC (validated routine, exact vs quantum) ----------
HBAR, KB, QE = 1.054571817e-34, 1.380649e-23, 1.602176634e-19
T = numpy.array(pts.evaluate())
E = numpy.array(pts.energies().inUnitsOf(eV))
Ej = E * QE
x = numpy.where(Ej / (KB * 300.0) < 1e-8, 1e-8, Ej / (KB * 300.0))
ex = numpy.exp(numpy.clip(x, 0, 700))
G = numpy.trapz(KB * T * (x ** 2) * ex / (ex - 1.0) ** 2, Ej / HBAR) / (2 * numpy.pi)
area = abs(numpy.cross(v1[:2], v2[:2])) * 1e-20
print("TRANS_MAX %.4f" % T.max())
print("G_total %.6e W/K   AREA %.4e m^2" % (G, area))
print("TBC %.4e W/m^2/K = %.2f MW/m^2/K" % (G / area, G / area / 1e6))
print("DEVTX_DONE %s" % TMD)
