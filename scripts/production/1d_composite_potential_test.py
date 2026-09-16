# PHASE 1d - CRITICAL FEASIBILITY TEST
# Can TremoloX build a COMPOSITE potential:
#     StillingerWeber_MoS_2013  (MoS2 intra)
#   + QuantumATK_MTP_HfO2_2022  (HfO2 intra, machine-learned)
#   + Lennard-Jones cross-terms (Mo-Hf, Mo-O, S-Hf, S-O)   <- our derived UFF/LB values
# and evaluate it on a stacked MoS2-on-HfO2 structure?
#
# This is the single decisive unknown for Phase 1d. Classical SW and a machine-learned
# MTP may not merge. Each stage is guarded so we learn exactly where it breaks.
#
# Cross-terms from 1d_derive_lj_crossterms.py (validated against D1 Table S1 to <1.2%):
#   Mo-Hf  eps=2.7535 meV  sigma=2.6059 A     [new]
#   Mo-O   eps=2.5136 meV  sigma=2.9186 A     [matches D1's published 2.51 / 2.920]
#   S-Hf   eps=6.0908 meV  sigma=3.0438 A     [new]
#   S-O    eps=5.5601 meV  sigma=3.3565 A     [matches D1's published 5.58 / 3.385]
from QuantumATK import *
import numpy

CHI = 1.0          # D1's interaction scaling; sweep 1..4 later
R_CUT = 10.0       # Angstrom

CROSS = [
    ("Mo", "Hf", 2.7535, 2.6059),
    ("Mo", "O",  2.5136, 2.9186),
    ("S",  "Hf", 6.0908, 3.0438),
    ("S",  "O",  5.5601, 3.3565),
]

print("COMPOSITE_TEST_START")

# ---------- stage 1: the two intra-material potentials ----------
try:
    sw = StillingerWeber_MoS_2013()
    print("S1_SW_OK   elements:", sw.supportedElements())
except Exception as ex:
    print("S1_SW_FAIL", type(ex).__name__, ex); raise SystemExit

try:
    mtp = QuantumATK_MTP_HfO2_2022()
    print("S2_MTP_OK  elements:", mtp.supportedElements())
except Exception as ex:
    print("S2_MTP_FAIL", type(ex).__name__, ex); raise SystemExit

# ---------- stage 2: merge them ----------
merged = None
for how in ["merge", "mergePotentials"]:
    try:
        cand = StillingerWeber_MoS_2013()
        getattr(cand, how)(QuantumATK_MTP_HfO2_2022())
        merged = cand
        print("S3_MERGE_OK via %s()  elements: %s" % (how, cand.supportedElements()))
        break
    except Exception as ex:
        print("S3_MERGE_FAIL via %s(): %s: %s" % (how, type(ex).__name__, ex))

if merged is None:
    print("S3_MERGE_ALL_FAILED -> composite potential not constructible this way")
    raise SystemExit

# ---------- stage 3: add the LJ cross-terms ----------
n_added = 0
for a, b, eps_meV, sig in CROSS:
    try:
        lj = LennardJonesPotential(
            particleType1=ParticleIdentifier(a, []),
            particleType2=ParticleIdentifier(b, []),
            r_cut=R_CUT * Angstrom,
            sigma=sig * Angstrom,
            epsilon=CHI * eps_meV * 1e-3 * eV,
        )
        merged.addPotential(lj)
        n_added += 1
        print("S4_LJ_OK   %s-%s eps=%.4f meV sigma=%.4f A" % (a, b, CHI * eps_meV, sig))
    except Exception as ex:
        print("S4_LJ_FAIL %s-%s -> %s: %s" % (a, b, type(ex).__name__, ex))
print("S4_LJ_ADDED %d of %d" % (n_added, len(CROSS)))

# ---------- stage 4: a stacked test structure ----------
# NOT physically optimised -- purely to check the potential EVALUATES on a cell that
# contains Mo, S, Hf and O simultaneously with a vdW gap.
try:
    a_hex = 3.18
    lat = UnitCell([a_hex * 3, 0.0, 0.0] * Angstrom,
                   [0.0, a_hex * 3, 0.0] * Angstrom,
                   [0.0, 0.0, 26.0] * Angstrom)
    els, xyz = [], []
    # a small HfO2 block (cubic fluorite fragment) near z = 3-8 Ang
    aH = 5.08
    for i in range(2):
        for j in range(2):
            x, y = i * aH, j * aH
            els += [Hafnium, Oxygen, Oxygen]
            xyz += [[x, y, 5.0],
                    [x + aH * 0.25, y + aH * 0.25, 5.0 + aH * 0.25],
                    [x + aH * 0.75, y + aH * 0.75, 5.0 - aH * 0.25]]
    # a MoS2 monolayer ~3 A above the HfO2 top surface
    z_mo = 5.0 + aH * 0.25 + 3.0 + 1.5645
    for i in range(3):
        for j in range(3):
            x = (i + 1.0 / 3.0) * a_hex
            y = (j + 2.0 / 3.0) * a_hex
            els += [Molybdenum, Sulfur, Sulfur]
            xyz += [[x, y, z_mo],
                    [x + a_hex / 3.0, y - a_hex / 3.0, z_mo + 1.5645],
                    [x + a_hex / 3.0, y - a_hex / 3.0, z_mo - 1.5645]]

    conf = BulkConfiguration(bravais_lattice=lat, elements=els,
                             cartesian_coordinates=numpy.array(xyz) * Angstrom)
    print("S5_STRUCT_OK  %d atoms (Mo/S/Hf/O)" % len(els))

    conf.setCalculator(TremoloXCalculator(parameters=merged))
    conf.update()
    e = TotalEnergy(conf).evaluate().inUnitsOf(eV)
    print("S6_EVAL_OK  E_tot = %.4f eV  (%.4f eV/atom)" % (e, e / len(els)))
    print("COMPOSITE_POTENTIAL_WORKS")
except Exception as ex:
    print("S5_S6_FAIL", type(ex).__name__, ex)

print("COMPOSITE_TEST_DONE")
