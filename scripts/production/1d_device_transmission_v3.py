"""PHASE 1d-C v2 - device phonon transmission -> TBC.

FIXES vs v1 (both real bugs, both caught by the sanity checks):
 1. INTERLAYER SPACING. v1 used d_tmd = C_TMD/2 = 3.075 A for MoS2. Bulk 2H-MoS2 has
    c = 12.29 A with TWO layers per cell, so the INTERLAYER SPACING IS 6.15 A. At 3.075 A
    the layers interpenetrate (a single layer is already 2*dz = 3.13 A thick) -> MIN_DISTANCE
    1.838 A and E = +1.09 eV/atom (positive = overlap).
 2. ELECTRODE CELL. QuantumATK requires the electrode's THIRD lattice vector to lie
    entirely along z. A cleaved (111) cell carries in-plane components in its c-vector
    (it is the ABC stacking vector). Both electrodes are therefore rebuilt explicitly with
    c = [0,0,L] and the ABC / 2H registry carried by the ATOM POSITIONS instead.

    fluorite (111): ABC stacking at in-plane fractions (0,0), (1/3,2/3), (2/3,1/3),
                    layer spacing d111 = a_cubic/sqrt(3); one electrode period = 3*d111
    2H TMD        : two layers per period, interlayer spacing c/2, registry alternates

usage: atkpython 1d_device_transmission_v2.py MoS2|WSe2 [n_ox_periods] [n_tmd_layers]
"""
import sys
import numpy
from QuantumATK import *
from SMW import *

TMD = sys.argv[1] if len(sys.argv) > 1 else 'MoS2'
N_OXP = int(sys.argv[2]) if len(sys.argv) > 2 else 2     # oxide ABC periods in central region
N_TMD = int(sys.argv[3]) if len(sys.argv) > 3 else 4     # TMD layers in central region

# a_TMD, dz, metal, chalcogen, bulk c (2 layers), relaxed vdW gap
_P = {'MoS2': (3.16, 1.5645, Molybdenum, Sulfur,   12.29, 4.697),
      'WSe2': (3.28, 1.67,   Tungsten,   Selenium, 12.98, 4.193)}
A_TMD, DZ, M_AT, X_AT, C_BULK, GAP = _P[TMD]
D_INTER = C_BULK / 2.0                                   # 6.145 / 6.49 A  <- the v1 bug
print("DEVTX2_START %s  interlayer=%.3f A  gap=%.3f A" % (TMD, D_INTER, GAP))

# ---------- oxide primitive (111) layer ----------
mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (cubic)'](
    amorphous_density=9.390 * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2})
ox1 = mat.pristineConfiguration().cleave(1, 1, 1, layers=1, outer_atom_index=1)
ov = numpy.array(ox1.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
a_ox = numpy.linalg.norm(ov[0])
D111 = ov[2][2]
oxp = numpy.array(ox1.cartesianCoordinates().inUnitsOf(Angstrom))
oxp[:, 2] -= oxp[:, 2].min()
oxe = list(ox1.elements())
print("OX_LAYER atoms=%d a=%.4f d111=%.4f" % (len(oxe), a_ox, D111))

# ---------- in-plane match ----------
best = None
for m in range(1, 14):
    for n in range(1, 16):
        mis = abs(n * A_TMD - m * a_ox) / (m * a_ox)
        if mis < 0.02 and (best is None or mis < best[0]):
            best = (mis, n, m)
mis, N_T, M_O = best
u1, u2 = ov[0], ov[1]                       # oxide surface vectors (in-plane)
v1, v2 = u1 * M_O, u2 * M_O                 # common supercell
w1, w2 = v1 / N_T, v2 / N_T                 # TMD vectors
print("MATCH %dx%s vs %dx HfO2 -> %.3f%%" % (N_T, TMD, M_O, 100 * mis))

ABC = [numpy.zeros(3), (u1 + 2 * u2) / 3.0, (2 * u1 + u2) / 3.0]


def add_oxide(els, cart, z0, nlayers):
    for k in range(nlayers):
        sh = ABC[k % 3] + numpy.array([0.0, 0.0, z0 + k * D111])
        for i in range(M_O):
            for j in range(M_O):
                s2 = sh + i * u1 + j * u2
                for e, q in zip(oxe, oxp):
                    els.append(e); cart.append(list(q + s2))
    return z0 + (nlayers - 1) * D111 + oxp[:, 2].max()


def add_tmd(els, cart, z_first_metal, nlayers):
    for k in range(nlayers):
        zc = z_first_metal + k * D_INTER
        off = (w1 + w2) / 3.0 if k % 2 else numpy.zeros(3)
        for i in range(N_T):
            for j in range(N_T):
                pm = (i + 1/3.) * w1 + (j + 2/3.) * w2 + off
                ps = (i + 2/3.) * w1 + (j + 1/3.) * w2 + off
                els += [M_AT, X_AT, X_AT]
                cart += [[pm[0], pm[1], zc], [ps[0], ps[1], zc + DZ], [ps[0], ps[1], zc - DZ]]
    return z_first_metal + (nlayers - 1) * D_INTER + DZ


# ---------- central region: oxide (incl. electrode copy) + gap + TMD (incl. copy) ----------
els, cart = [], []
n_ox_layers = 3 * (N_OXP + 1)                      # +1 period -> electrode overlap
z_ox_top = add_oxide(els, cart, 0.0, n_ox_layers)
z_first_metal = z_ox_top + GAP + DZ
z_tmd_top = add_tmd(els, cart, z_first_metal, N_TMD + 2)
C = z_tmd_top + (D_INTER - DZ)                     # close the stack periodically
cell = UnitCell(v1 * Angstrom, v2 * Angstrom, [0.0, 0.0, C] * Angstrom)
central = BulkConfiguration(bravais_lattice=cell, elements=els,
                            cartesian_coordinates=numpy.array(cart) * Angstrom)
print("CENTRAL atoms=%d C=%.3f z_ox_top=%.3f" % (len(els), C, z_ox_top))

p = numpy.array(central.cartesianCoordinates().inUnitsOf(Angstrom))
dmin = 1e9
for i in range(0, len(p), 3):
    dv = p[i] - p
    dv[:, 2] -= C * numpy.round(dv[:, 2] / C)
    r = numpy.sqrt((dv ** 2).sum(axis=1)); r[r < 1e-6] = 1e9
    dmin = min(dmin, r.min())
print("MIN_DISTANCE %.3f A" % dmin)
if dmin < 1.2:
    print("SANITY_ABORT"); raise SystemExit

pot = TorchX_MACE_MPA_0_medium(); pot.merge(DispersionD3BJ())
calc = TremoloXCalculator(parameters=pot)
central.setCalculator(calc); central.update()
e = TotalEnergy(central).evaluate().inUnitsOf(eV)
print("E_central %.4f eV (%.4f eV/atom)" % (e, e / len(els)))
if e > 0:
    print("POSITIVE_ENERGY_ABORT"); raise SystemExit

# ---------- electrodes: EXTRACTED from the central region ----------
# v2 failed with "right electrode extension does not match the repeated right electrode":
# QuantumATK requires the central region to BEGIN and END with exact copies of the
# electrodes. Building them independently gave the right periodicity but the wrong PHASE
# (electrode metal 4.58 A below its cell top vs 6.15 A in the central region). Extracting
# the terminal slabs from the central region makes the match exact by construction.
allp = numpy.array(central.cartesianCoordinates().inUnitsOf(Angstrom))
alle = list(central.elements())

L_OX = 3 * D111          # one fluorite-(111) ABC period
L_TM = C_BULK            # one 2H period

# LEFT electrode = bottom L_OX of the central region
sel = [i for i in range(len(alle)) if allp[i][2] < L_OX - 1e-6]
oe = [alle[i] for i in sel]
oc = [[allp[i][0], allp[i][1], allp[i][2]] for i in sel]
ox_elec = BulkConfiguration(
    bravais_lattice=UnitCell(v1 * Angstrom, v2 * Angstrom, [0.0, 0.0, L_OX] * Angstrom),
    elements=oe, cartesian_coordinates=numpy.array(oc) * Angstrom)

# RIGHT electrode = top L_TM of the central region, shifted to start at z=0
sel = [i for i in range(len(alle)) if allp[i][2] > C - L_TM - 1e-6]
te = [alle[i] for i in sel]
tc = [[allp[i][0], allp[i][1], allp[i][2] - (C - L_TM)] for i in sel]
tmd_elec = BulkConfiguration(
    bravais_lattice=UnitCell(v1 * Angstrom, v2 * Angstrom, [0.0, 0.0, L_TM] * Angstrom),
    elements=te, cartesian_coordinates=numpy.array(tc) * Angstrom)

print("ELECTRODES oxide=%d (c=%.3f)   %s=%d (c=%.3f)  [extracted from central]"
      % (len(oe), L_OX, TMD, len(te), L_TM))

dev = DeviceConfiguration(central, [ox_elec, tmd_elec])
dev.setCalculator(calc); dev.update()
print("DEVICE_OK")

dm = DynamicalMatrix(configuration=dev, filename='tx2_%s.hdf5' % TMD, object_id='dm',
                     calculator=calc, repetitions=(1, 1, 1),
                     use_equivalent_bulk=False, constrain_electrodes=False)
dm.update(); print("DM_OK")

pts = PhononTransmissionSpectrum(
    dynamical_matrix=dm, energies=numpy.linspace(0, 0.08, 161) * eV,
    qpoints=MonkhorstPackGrid(5, 5), infinitesimal=1e-6 * eV,
    self_energy_calculator=RecursionSelfEnergy())
nlsave('tx2_%s.hdf5' % TMD, pts); print("PTS_OK")

HBAR, KB, QE = 1.054571817e-34, 1.380649e-23, 1.602176634e-19
Tr = numpy.array(pts.evaluate()); E = numpy.array(pts.energies().inUnitsOf(eV))
Ej = E * QE
x = numpy.where(Ej / (KB * 300.0) < 1e-8, 1e-8, Ej / (KB * 300.0))
ex = numpy.exp(numpy.clip(x, 0, 700))
G = numpy.trapz(KB * Tr * (x ** 2) * ex / (ex - 1.0) ** 2, Ej / HBAR) / (2 * numpy.pi)
area = abs(v1[0] * v2[1] - v1[1] * v2[0]) * 1e-20
print("TRANS_MAX %.4f" % Tr.max())
print("G_total %.6e W/K  AREA %.4e m2" % (G, area))
print("TBC %.2f MW/m^2/K   (D1 MoS2/SiO2 ref: 11.5 chi=1, 25.6 chi=2, 13.0 saturated)"
      % (G / area / 1e6))
print("DEVTX2_DONE %s" % TMD)
