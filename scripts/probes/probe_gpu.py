from QuantumATK import *
import inspect
pot = TorchX_MACE_MPA_0_medium(); pot.merge(DispersionD3BJ())
calc = TremoloXCalculator(parameters=pot)
for m in ('getGPUDevices','setGPUDevices','setGPUDevice','setGPUAcceleration','setAllGPUAcceleration'):
    f = getattr(calc, m, None)
    try:
        print("%-22s %s" % (m, inspect.signature(f)))
    except Exception as e:
        print("%-22s (sig err %s)" % (m, e))
    d = inspect.getdoc(f)
    if d:
        print("     doc:", d.strip().splitlines()[0][:150])
print("--- current ---")
try:
    print("getGPUDevices() ->", calc.getGPUDevices())
except Exception as e:
    print("getGPUDevices err", e)
for m in ('_gpuEnabledAndUsed','_gpuTorchXEnabled'):
    try:
        print(m, "->", getattr(calc, m)())
    except Exception as e:
        print(m, "err", e)
print("--- try both GPUs ---")
try:
    calc.setGPUDevices([0, 1]); print("setGPUDevices([0,1]) OK ->", calc.getGPUDevices())
except Exception as e:
    print("setGPUDevices([0,1]) FAILED:", type(e).__name__, e)
print("PROBE2_DONE")
