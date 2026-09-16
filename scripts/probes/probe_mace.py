from QuantumATK import *
import inspect
pot = TorchX_MACE_MPA_0_medium()
print("POT_TYPE", type(pot).__name__)
try:
    print("SIG", inspect.signature(TorchX_MACE_MPA_0_medium))
except Exception as e:
    print("SIG_ERR", e)
print("--- attrs mentioning dtype/precision/device/batch/chunk ---")
for a in dir(pot):
    if any(k in a.lower() for k in ('dtype','float','precis','device','batch','chunk','cuda','gpu','memory','double','single')):
        print("  POT.", a)
print("--- TremoloXCalculator signature ---")
try:
    print(inspect.signature(TremoloXCalculator))
except Exception as e:
    print("ERR", e)
for a in dir(TremoloXCalculator):
    if any(k in a.lower() for k in ('dtype','float','precis','device','batch','chunk','cuda','gpu','memory')):
        print("  CALC.", a)
print("--- other MACE/TorchX potentials available ---")
import QuantumATK as Q
names = [n for n in dir(Q) if 'MACE' in n.upper() or n.startswith('TorchX')]
for n in names:
    print("  ", n)
print("PROBE_DONE")
