"""PHASE 1d-B (v10) - MoS2 / crystalline HfO2(111), following G3.

WHY CRYSTALLINE (after 8 failed amorphous builds):
  G3 = KC, Longo, Wallace, Cho, ACS Omega 2, 2827 (2017) explicitly justifies a
  CRYSTALLINE HfO2 slab as the model for an ALD oxide:
    "the local Hf-O bonding is more important than long-range order for interface
     engineering. Therefore, the 15 dielectric layers considered represent a good
     model system for this type of interfaces."
  The amorphous requirement was MY inference from E4, not a constraint the field imposes.

G3 RECIPE (copied):
  * S-terminated MoS2 surface on an O-terminated HfO2 (111) surface
  * lattice mismatch < 1%   (ours: 8xMoS2 / 7xHfO2(111) = 0.54%)
  * 5 atomic layers of Hf + 10 of O
  * 16 A vacuum
  * weak vdW coupling, NO covalent bonding (their Fig. 1c: flat vacuum level)

CRITICAL: build the surface with BulkConfiguration.cleave(), NOT by hand.
v1 failed precisely because I hand-placed the fluorite O sublattice (both O atoms on the
Hf site -> 0.733 A contacts -> +39.5 eV/atom).

usage: atkpython 1d_iface_v10_crystalline.py [n_layers] [gap] [vacuum]
"""
import sys
import numpy
from QuantumATK import *
from SMW import *

NLAYERS = int(sys.argv[1]) if len(sys.argv) > 1 else 5      # Hf layers (G3: 5)
GAP = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
VACUUM = float(sys.argv[3]) if len(sys.argv) > 3 else 16.0  # G3: 16 A

A_MOS2, DZ_S = 3.16, 1.5645

print("V10_START  crystalline HfO2(111) route, G3 recipe")

# ---------- 1. cubic HfO2 from the vendor database, then CLEAVE (111) ----------
mat = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (cubic)'](
    amorphous_density=9.390 * gram / cm ** 3,
    elements_with_stoichiometry={Hafnium: 1, Oxygen: 2})
bulk = mat.pristineConfiguration()
cb = numpy.array(bulk.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
print("CUBIC_HFO2 atoms=%d  a=%.4f b=%.4f c=%.4f"
      % (len(bulk.elements()),
         numpy.linalg.norm(cb[0]), numpy.linalg.norm(cb[1]), numpy.linalg.norm(cb[2])))

surf = None
# outer_atom_index=1 gives an O-TERMINATED surface on BOTH faces (probe_termination.py).
# v10 used the default (oai=0) -> Hf-terminated. A bare Hf surface is reactive: it grabbed
# the MoS2 sulfur and collapsed the vdW gap from 3.0 to 0.81 A, distorting the monolayer
# to 5.0 A spread. G3 specifies O-termination precisely because it is passivated and gives
# the weak vdW interaction (their flat vacuum level, no charge transfer).
for kwargs in [dict(layers=NLAYERS, outer_atom_index=1,
                    top_vacuum=VACUUM * Angstrom, bottom_vacuum=2.0 * Angstrom),
               dict(layers=NLAYERS, outer_atom_index=1)]:
    try:
        surf = bulk.cleave(1, 1, 1, **kwargs)
        print("CLEAVE_OK with %s" % list(kwargs.keys()))
        break
    except Exception as ex:
        print("  cleave(%s) failed: %s: %s" % (list(kwargs.keys()), type(ex).__name__, str(ex)[:150]))
if surf is None:
    print("CLEAVE_ALL_FAILED")
    raise SystemExit

sv = numpy.array(surf.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
a_s, b_s = numpy.linalg.norm(sv[0]), numpy.linalg.norm(sv[1])
ang = numpy.degrees(numpy.arccos(numpy.dot(sv[0], sv[1]) / (a_s * b_s)))
sp = numpy.array(surf.cartesianCoordinates().inUnitsOf(Angstrom))
ss = [str(e.symbol()) for e in surf.elements()]
print("SURFACE  atoms=%d  a=%.4f b=%.4f  angle=%.2f deg  c=%.3f"
      % (len(ss), a_s, b_s, ang, numpy.linalg.norm(sv[2])))
print("SURFACE_Z  min %.3f max %.3f   (Hf %d, O %d)"
      % (sp[:, 2].min(), sp[:, 2].max(), ss.count("Hf"), ss.count("O")))
# termination check: what is the topmost species?
ztop = sp[:, 2].max()
top_sp = [ss[i] for i in range(len(ss)) if sp[i, 2] > ztop - 0.6]
print("TOP_TERMINATION %s  (G3 wants O-terminated)" % sorted(set(top_sp)))

# ---------- 2. lattice match: n*a_MoS2 ~ m*a_surf ----------
best = None
for m in range(1, 13):
    for n in range(1, 15):
        mis = abs(n * A_MOS2 - m * a_s) / (m * a_s)
        if mis < 0.02 and (best is None or mis < best[0]):
            best = (mis, n, m)
mis, N_MOS2, M_HF = best
print("LATTICE_MATCH  %d x MoS2 (%.3f) vs %d x HfO2surf (%.3f)  ->  mismatch %.3f%%"
      % (N_MOS2, N_MOS2 * A_MOS2, M_HF, M_HF * a_s, 100 * mis))
if mis > 0.01:
    print("MISMATCH_TOO_BIG (G3 requires <1%)")

slab = surf.repeat(M_HF, M_HF, 1)
sp = numpy.array(slab.cartesianCoordinates().inUnitsOf(Angstrom))
els = list(slab.elements())
lv = numpy.array(slab.bravaisLattice().primitiveVectors().inUnitsOf(Angstrom))
L1, L2 = lv[0], lv[1]
print("SLAB_REPEATED atoms=%d" % len(els))

# ---------- 3. put MoS2 on top, strained onto the HfO2 surface lattice ----------
z_top = sp[:, 2].max()
z_mo = z_top + GAP + DZ_S
# MoS2 lattice vectors are the SLAB vectors divided by N_MOS2 (so MoS2 conforms exactly)
w1 = L1 / N_MOS2
w2 = L2 / N_MOS2
cart = [list(q) for q in sp]
for i in range(N_MOS2):
    for j in range(N_MOS2):
        pm = (i + 1.0 / 3.0) * w1 + (j + 2.0 / 3.0) * w2
        ps = (i + 2.0 / 3.0) * w1 + (j + 1.0 / 3.0) * w2
        els += [Molybdenum, Sulfur, Sulfur]
        cart += [[pm[0], pm[1], z_mo],
                 [ps[0], ps[1], z_mo + DZ_S],
                 [ps[0], ps[1], z_mo - DZ_S]]

C = z_mo + DZ_S + VACUUM
cell = UnitCell(L1 * Angstrom, L2 * Angstrom, [0.0, 0.0, C] * Angstrom)
iface = BulkConfiguration(bravais_lattice=cell, elements=els,
                          cartesian_coordinates=numpy.array(cart) * Angstrom)
n_tot = len(els)
print("IFACE_ATOMS %d  (HfO2 %d + MoS2 %d)  C=%.2f  z_top=%.2f  z_Mo=%.2f"
      % (n_tot, len(sp), n_tot - len(sp), C, z_top, z_mo))

# ---------- 4. sanity checks BEFORE any relaxation ----------
p = numpy.array(iface.cartesianCoordinates().inUnitsOf(Angstrom))
sy = [str(e.symbol()) for e in iface.elements()]
dmin = 1e9
for i in range(0, len(p), 2):
    dv = p[i] - p
    dv[:, 0] -= L1[0] * numpy.round(dv[:, 0] / L1[0])
    r = numpy.sqrt((dv ** 2).sum(axis=1))
    r[r < 1e-6] = 1e9
    dmin = min(dmin, r.min())
print("MIN_DISTANCE %.3f A" % dmin)
# Mo-S bond check
mo = [i for i, t in enumerate(sy) if t == "Mo"]
su = [i for i, t in enumerate(sy) if t == "S"]
bb = 1e9
for i in mo[:6]:
    dv = p[i] - p[su]
    r = numpy.sqrt((dv ** 2).sum(axis=1))
    bb = min(bb, r.min())
print("MOS2_BOND %.3f A  (must be ~2.41)" % bb)
if dmin < 1.0 or abs(bb - 2.41) > 0.12:
    print("SANITY_ABORT")
    raise SystemExit

# ---------- 5. relax ----------
pot = TorchX_MACE_MPA_0_medium()
pot.merge(DispersionD3BJ())
iface.setCalculator(TremoloXCalculator(parameters=pot))
iface.update()
e0 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_BEFORE %.4f eV (%.4f eV/atom)" % (e0, e0 / n_tot))

iface = OptimizeGeometry(iface, max_forces=0.02 * eV / Ang, max_steps=400,
                         optimize_cell=False, trajectory_filename='v10_relax.hdf5')
e1 = TotalEnergy(iface).evaluate().inUnitsOf(eV)
print("E_AFTER  %.4f eV (%.4f eV/atom)" % (e1, e1 / n_tot))

q = numpy.array(iface.cartesianCoordinates().inUnitsOf(Angstrom))
z_ox = max(q[i][2] for i, t in enumerate(sy) if t in ("Hf", "O"))
z_s = min(q[i][2] for i, t in enumerate(sy) if t == "S")
mz = [q[i][2] for i, t in enumerate(sy) if t in ("Mo", "S")]
print("FINAL_GAP   %.3f A  (target ~3, G3 vdW interface)" % (z_s - z_ox))
print("MOS2_SPREAD %.3f A  (ideal %.3f)" % (max(mz) - min(mz), 2 * DZ_S))
nlsave('mos2_hfo2_v10.hdf5', iface)
print("SAVED mos2_hfo2_v10.hdf5")
print("V10_DONE")
