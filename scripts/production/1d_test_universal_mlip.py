# PHASE 1d - can a UNIVERSAL machine-learned potential replace the SW+MTP+LJ stack?
#
# Discovery (2026-09-03): QuantumATK Y-2026.03 ships TorchX foundation MLIPs trained on
# Materials Project that support ~89 elements INCLUDING Mo, S, W, Se, Hf, O.
# If one of these reproduces our DFT phonons it would:
#   * solve the WSe2 force-field BLOCKER outright,
#   * remove the need for hand-derived LJ cross-terms (vdW handled natively),
#   * give ONE consistent potential for MoS2, WSe2, HfO2 and every interface between them.
#
# GROUND TRUTH from our own DFT (Phase 1b/1c), same geometry, same settings:
#   MoS2 : max 457.4 cm-1 ; Gamma 276.9(x2) 373.6(x2) 393.2 457.4
#   WSe2 : max 298.2 cm-1 ; Gamma 166.1 166.3 231.9 237.6 237.8 298.2
# Literature: MoS2 max ~470-477 ; WSe2 max ~308
from QuantumATK import *
import numpy

CM1 = 8065.54
MATS = {
    "MoS2": (3.18, 1.5645, [Molybdenum, Sulfur, Sulfur], 457.4),
    "WSe2": (3.28, 1.67,   [Tungsten, Selenium, Selenium], 298.2),
}
CANDIDATES = [
    "TorchX_MACE_MPA_0_medium",
    "TorchX_MACE_MP_0_L1_2023",
    "TorchX_M3GNet_MP_PES_2021",
    "TorchX_CHGNet_0_3_0",
    "TorchX_MatterSim_v1_0_0_5M",
]


def build(mat):
    a, dz, els, _ = MATS[mat]
    c = 20.0
    f = dz / c
    return BulkConfiguration(
        bravais_lattice=Hexagonal(a * Angstrom, c * Angstrom),
        elements=els,
        fractional_coordinates=[[1/3., 2/3., 0.50],
                                [2/3., 1/3., 0.50 + f],
                                [2/3., 1/3., 0.50 - f]])


print("MLIP_TEST_START")

working = []
for pname in CANDIDATES:
    print("=" * 72)
    print("CANDIDATE", pname)
    # stage 1: does it even instantiate + evaluate? (TorchX models may need a venv)
    try:
        pot = globals()[pname]()
        calc = TremoloXCalculator(parameters=pot)
        conf = build("MoS2")
        conf.setCalculator(calc)
        conf.update()
        e = TotalEnergy(conf).evaluate().inUnitsOf(eV)
        print("  EVAL_OK  E(MoS2) = %.4f eV" % e)
    except Exception as ex:
        print("  EVAL_FAIL %s: %s" % (type(ex).__name__, str(ex)[:220]))
        continue

    # stage 2: phonons for both materials vs our DFT ground truth
    ok = True
    for mat in ["MoS2", "WSe2"]:
        try:
            conf = build(mat)
            conf.setCalculator(TremoloXCalculator(parameters=globals()[pname]()))
            conf.update()
            conf = OptimizeGeometry(conf, max_forces=0.02*eV/Ang, max_steps=100,
                                    optimize_cell=False)
            dm = DynamicalMatrix(configuration=conf,
                                 filename='mlip_%s_%s.hdf5' % (pname, mat),
                                 object_id='dm',
                                 calculator=TremoloXCalculator(parameters=globals()[pname]()),
                                 repetitions=(3, 3, 1))
            dm.update()
            pb = PhononBandstructure(dynamical_matrix=dm, route=['G', 'M', 'K', 'G'],
                                     points_per_segment=30)
            w = numpy.array(pb.evaluate().inUnitsOf(eV)) * CM1
            ref = MATS[mat][3]
            g = numpy.sort(w[0, :])
            print("  %s  max=%7.1f  (our DFT %6.1f, err %+6.1f%%)  acoustic=%s"
                  % (mat, w.max(), ref, 100*(w.max()-ref)/ref,
                     " ".join("%.1f" % x for x in g[:3])))
        except Exception as ex:
            print("  %s PHONON_FAIL %s: %s" % (mat, type(ex).__name__, str(ex)[:180]))
            ok = False
    if ok:
        working.append(pname)

print("=" * 72)
print("WORKING_MLIPS", working)
print("MLIP_TEST_DONE")
