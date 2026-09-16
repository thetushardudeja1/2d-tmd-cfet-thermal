# Probe: what force fields / potentials ship with QuantumATK for TMD + HfO2 elements,
# and confirm the phonon API objects import.
import re

print("PROBE_START")

# 1) force-field library
try:
    import AddOns.ForceFields as FF
    names = [n for n in dir(FF) if not n.startswith("_")]
    print("FF_TOTAL", len(names))
    pat = re.compile(r"MoS|WSe|WS2|MoSe|Hf|Stillinger|Tersoff|SW|TMD|Chalco|Brenner|REBO", re.I)
    hits = [n for n in names if pat.search(n)]
    print("FF_HITS", len(hits))
    for h in hits[:40]:
        print("   FF:", h)
except Exception as ex:
    print("FF_IMPORT_FAILED", type(ex).__name__, ex)

# 2) phonon / transport API availability
from QuantumATK import *
for obj in ["DynamicalMatrix", "PhononBandstructure", "PhononDensityOfStates",
            "PhononTransmissionSpectrum", "ThermalConductivity", "MTPPotential",
            "TremoloXCalculator", "LCAOCalculator"]:
    print("API %-28s %s" % (obj, "OK" if obj in globals() else "MISSING"))

print("PROBE_OK")
