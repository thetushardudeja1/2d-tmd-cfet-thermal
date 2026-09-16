"""PHASE 1d-B (v2) - generate AMORPHOUS HfO2 by melt-quench.

WHY AMORPHOUS (not the crystalline slab of v1):
  * Literature: E4 (npj 2D Mater. 2025) - real MoS2 gate stacks use **amorphous ALD HfO2**.
    25 ALD mentions, explicit "amorphous HfO2", zero monoclinic/polycrystalline.
  * It removes the lattice-matching problem entirely (no commensuration, no strain).
  * D1 used amorphous SiO2 under MoS2 - same precedent.
  * v1 hand-built fluorite(111) was geometrically WRONG (see INTERFACE_POTENTIAL_RECIPE.md).

TEMPLATE: vendor SMW example `HKMG_builder_example.py` (High-K Metal Gate).
CALCULATOR: QuantumATK_MTP_HfO2_2022 - purpose-trained for HfO2 and far cheaper than MACE
            for a 40 ps melt-quench. The interface relaxation later uses MACE_MPA + D3BJ.

usage: atkpython 1d_build_amorphous_hfo2.py [length_A] [melt_ps]
"""
import sys
from QuantumATK import *
from SMW import *

LENGTH = float(sys.argv[1]) if len(sys.argv) > 1 else 15.0    # slab thickness, Angstrom
MELT_PS = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0   # melt time, ps

print("AMORPH_START length=%.1f A melt=%.1f ps" % (LENGTH, MELT_PS))

# cheap, purpose-trained potential for the melt-quench
mtp_calc = TremoloXCalculator(parameters=QuantumATK_MTP_HfO2_2022())

amorphous_method = MolecularDynamicsMeltQuench(
    melt_temperature=3500 * Kelvin,
    melt_time=MELT_PS * ps,
    quench_temperature=300 * Kelvin,
    quench_rate=74 * Kelvin / ps,
    time_step=1.0 * fs,
    max_retries=30,
    random_seed=1234,
)

try:
    material = MaterialSpecificationsDatabase.HKMG_MATERIALS['Hafnia (monoclinic)'](
        amorphous_density=9.390 * gram / cm ** 3,
        elements_with_stoichiometry={Hafnium: 1, Oxygen: 2},
        amorphization_calculator=mtp_calc,
    )
    print("MATERIAL_OK  Hafnia, rho=9.390 g/cm3, Hf:O = 1:2")
except Exception as ex:
    print("MATERIAL_FAIL %s: %s" % (type(ex).__name__, str(ex)[:300]))
    raise SystemExit

try:
    layer = AmorphousLayer(material, length=LENGTH * Angstrom, method=amorphous_method)
    print("LAYER_OK")
except Exception as ex:
    print("LAYER_FAIL %s: %s" % (type(ex).__name__, str(ex)[:300]))
    raise SystemExit

# Realise the amorphous structure. API surface varies; try the plausible entry points
# and report which one worked so the recipe can be recorded.
conf = None
for meth in ["configuration", "buildConfiguration", "build", "update", "generate"]:
    if not hasattr(layer, meth):
        continue
    try:
        out = getattr(layer, meth)()
        if out is not None:
            conf = out
            print("REALISED_VIA %s()" % meth)
            break
        else:
            print("CALLED %s() -> None" % meth)
    except Exception as ex:
        print("  %s() failed: %s: %s" % (meth, type(ex).__name__, str(ex)[:160]))

if conf is None:
    print("NO_CONFIG - available methods:", [m for m in dir(layer) if not m.startswith("_")][:30])
    raise SystemExit

try:
    n = len(conf.elements())
    print("AMORPH_ATOMS %d" % n)
    conf.setCalculator(mtp_calc)
    conf.update()
    e = TotalEnergy(conf).evaluate().inUnitsOf(eV)
    print("AMORPH_E %.4f eV (%.4f eV/atom)" % (e, e / n))
    nlsave('amorphous_hfo2.hdf5', conf)
    print("SAVED amorphous_hfo2.hdf5")
except Exception as ex:
    print("EVAL_FAIL %s: %s" % (type(ex).__name__, str(ex)[:250]))

print("AMORPH_DONE")
