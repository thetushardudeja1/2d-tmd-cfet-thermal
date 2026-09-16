# Last GPU path: the API-level gpu_acceleration flag (vendor uses gpu_acceleration=Automatic
# in training_example.py). Does it work under PLAIN atkpython, without atkpython_gpu?
from QuantumATK import *

print("GPUAPI_START")
for name in ["Enabled", "Disabled", "Automatic"]:
    print("SYMBOL %-10s %s" % (name, "present" if name in globals() else "MISSING"))

# ParallelConjugateGradientSolver documents gpu_acceleration=Enabled|Disabled|Automatic
for mode in ["Automatic", "Enabled"]:
    try:
        s = ParallelConjugateGradientSolver(gpu_acceleration=globals()[mode])
        print("SOLVER_OK  gpu_acceleration=%s -> constructed" % mode)
    except Exception as ex:
        print("SOLVER_FAIL gpu_acceleration=%s -> %s: %s" % (mode, type(ex).__name__, ex))

# Does a real calculation using that solver run, and does it touch the GPU?
try:
    lattice = FaceCenteredCubic(5.4306 * Angstrom)
    conf = BulkConfiguration(bravais_lattice=lattice, elements=[Silicon, Silicon],
                             fractional_coordinates=[[0., 0., 0.], [.25, .25, .25]])
    nap = NumericalAccuracyParameters(k_point_sampling=(3, 3, 3))
    conf.setCalculator(LCAOCalculator(numerical_accuracy_parameters=nap))
    e = TotalEnergy(conf).evaluate().inUnitsOf(eV)
    print("SCF_RAN_OK E=%.3f eV" % e)
except Exception as ex:
    print("SCF_FAIL %s: %s" % (type(ex).__name__, ex))

print("GPUAPI_OK")
