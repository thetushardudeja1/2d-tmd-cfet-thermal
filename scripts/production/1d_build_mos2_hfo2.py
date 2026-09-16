"""PHASE 1d-B - build and relax the MoS2 / HfO2 interface.

LATTICE MATCH (reported honestly in the paper):
  MoS2 monolayer, hexagonal          a = 3.16 A
  cubic fluorite HfO2 (111) surface  a_surf = a_cubic/sqrt(2) = 5.08/1.4142 = 3.5921 A
  8 x MoS2 = 25.280 A ; 7 x HfO2(111) = 25.145 A  ->  MISMATCH 0.54%
  We strain HfO2 by +0.54% onto the MoS2 lattice. (Alternatives were far worse:
  6/5 -> 5.6%, 5/4 -> 9.9%.)

FLUORITE (111) GEOMETRY (derived, not guessed):
  Hf planes repeat along [111] with spacing a/sqrt(3) = 2.9330 A
  O planes sit at +-a*sqrt(3)/12 = +-0.7333 A about each Hf plane (O-Hf-O trilayer)
  in-plane ABC stacking at fractional (0,0), (1/3,2/3), (2/3,1/3)

POTENTIAL: MACE_MPA + D3BJ (validated: phonons <1%, vdW binding -24.94 meV/atom @ 6.40 A)

usage: atkpython 1d_build_mos2_hfo2.py [n_hf_layers] [gap_A]
"""
import sys
import numpy
from QuantumATK import *

N_HF = int(sys.argv[1]) if len(sys.argv) > 1 else 3      # Hf planes in the slab
GAP = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0   # initial vdW gap (A)

A_MOS2 = 3.16
DZ_S = 1.5645
N_MOS2 = 8            # MoS2 supercell
A_HFO2_CUBIC = 5.08
N_HF_CELL = 7         # HfO2(111) supercell

A_SURF = A_HFO2_CUBIC / numpy.sqrt(2.0)          # 3.5921
D111 = A_HFO2_CUBIC / numpy.sqrt(3.0)            # 2.9330  Hf-Hf along [111]
DO = A_HFO2_CUBIC * numpy.sqrt(3.0) / 12.0       # 0.7333  O offset about Hf

L = N_MOS2 * A_MOS2                              # common in-plane length 25.280
strain = (L / (N_HF_CELL * A_SURF) - 1.0) * 100.0
a_hf_strained = L / N_HF_CELL

print("BUILD_START")
print("IN_PLANE_L      %.4f A" % L)
print("HFO2_STRAIN     %+.3f %%  (a_surf %.4f -> %.4f A)" % (strain, A_SURF, a_hf_strained))

C = 60.0        # tall cell: slab + gap + MoS2 + vacuum
els, cart = [], []

# ---------------- HfO2(111) slab, ABC stacked, O-Hf-O trilayers ----------------
# hexagonal in-plane basis for the HfO2 sublattice (strained to a_hf_strained)
ah = a_hf_strained
v1 = numpy.array([ah, 0.0])
v2 = numpy.array([-ah / 2.0, ah * numpy.sqrt(3.0) / 2.0])
abc = [numpy.array([0.0, 0.0]),
       numpy.array([1.0 / 3.0, 2.0 / 3.0]),
       numpy.array([2.0 / 3.0, 1.0 / 3.0])]

z0 = 8.0
for k in range(N_HF):
    zk = z0 + k * D111
    off = abc[k % 3]
    for i in range(N_HF_CELL):
        for j in range(N_HF_CELL):
            p = (i + off[0]) * v1 + (j + off[1]) * v2
            els.append(Hafnium); cart.append([p[0], p[1], zk])
            els.append(Oxygen);  cart.append([p[0], p[1], zk + DO])
            els.append(Oxygen);  cart.append([p[0], p[1], zk - DO])
n_hfo2 = len(els)
z_slab_top = z0 + (N_HF - 1) * D111 + DO
print("HFO2_ATOMS      %d  (%d Hf planes, slab top z = %.3f A)" % (n_hfo2, N_HF, z_slab_top))

# ---------------- MoS2 monolayer above the slab ----------------
am = A_MOS2
w1 = numpy.array([am, 0.0])
w2 = numpy.array([-am / 2.0, am * numpy.sqrt(3.0) / 2.0])
z_mo = z_slab_top + GAP + DZ_S          # gap measured slab-top -> lower S plane
for i in range(N_MOS2):
    for j in range(N_MOS2):
        pm = (i + 1.0 / 3.0) * w1 + (j + 2.0 / 3.0) * w2
        ps = (i + 2.0 / 3.0) * w1 + (j + 1.0 / 3.0) * w2
        els.append(Molybdenum); cart.append([pm[0], pm[1], z_mo])
        els.append(Sulfur);     cart.append([ps[0], ps[1], z_mo + DZ_S])
        els.append(Sulfur);     cart.append([ps[0], ps[1], z_mo - DZ_S])
n_mos2 = len(els) - n_hfo2
print("MOS2_ATOMS      %d" % n_mos2)
print("TOTAL_ATOMS     %d" % len(els))

lattice = UnitCell([L, 0.0, 0.0] * Angstrom,
                   [-L / 2.0, L * numpy.sqrt(3.0) / 2.0, 0.0] * Angstrom,
                   [0.0, 0.0, C] * Angstrom)
conf = BulkConfiguration(bravais_lattice=lattice, elements=els,
                         cartesian_coordinates=numpy.array(cart) * Angstrom)

# ---------------- potential: the validated production choice ----------------
pot = TorchX_MACE_MPA_0_medium()
pot.merge(DispersionD3BJ())
calc = TremoloXCalculator(parameters=pot)
conf.setCalculator(calc)
conf.update()
e0 = TotalEnergy(conf).evaluate().inUnitsOf(eV)
print("E_INITIAL       %.4f eV  (%.4f eV/atom)" % (e0, e0 / len(els)))

# ---------------- relax (D1 protocol: minimise forces) ----------------
try:
    conf = OptimizeGeometry(conf, max_forces=0.02 * eV / Ang, max_steps=300,
                            optimize_cell=False,
                            trajectory_filename='mos2_hfo2_relax.hdf5')
    e1 = TotalEnergy(conf).evaluate().inUnitsOf(eV)
    print("E_RELAXED       %.4f eV  (%.4f eV/atom)" % (e1, e1 / len(els)))
    print("E_DROP          %.4f eV" % (e1 - e0))
except Exception as ex:
    print("RELAX_FAIL %s: %s" % (type(ex).__name__, str(ex)[:200]))

# ---------------- report the final geometry ----------------
try:
    pos = conf.cartesianCoordinates().inUnitsOf(Angstrom)
    sym = [str(e.symbol()) for e in conf.elements()]
    z = numpy.array([p[2] for p in pos])
    z_hf_top = max(z[i] for i, s in enumerate(sym) if s in ("Hf", "O"))
    z_s_bot = min(z[i] for i, s in enumerate(sym) if s == "S")
    print("FINAL_GAP       %.3f A   (slab top -> lowest S)" % (z_s_bot - z_hf_top))
    print("MOS2_Z_SPREAD   %.3f A" % (max(z[i] for i, s in enumerate(sym) if s in ("Mo", "S"))
                                      - min(z[i] for i, s in enumerate(sym) if s in ("Mo", "S"))))
except Exception as ex:
    print("GEOM_REPORT_FAIL %s" % ex)

nlsave('mos2_hfo2_interface.hdf5', conf)
print("SAVED mos2_hfo2_interface.hdf5")
print("BUILD_DONE")
