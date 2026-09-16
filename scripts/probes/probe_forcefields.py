# Probe: which shipped TremoloX potentials cover the TMD / HfO2 element set.
# Vendor pattern (phonon_transmission.py):  potentialSet = StillingerWeber_Si_1985()
from QuantumATK import *
import re

g = sorted(globals().keys())
print("GLOBALS_TOTAL", len(g))

# potentials look like StillingerWeber_Si_1985, Tersoff_C_1989, etc.
potpat = re.compile(r"StillingerWeber|Tersoff|Brenner|REBO|Ternary|_\d{4}$|Potential$", re.I)
pots = [n for n in g if potpat.search(n)]
print("POTENTIAL_LIKE_COUNT", len(pots))

# now filter for our elements
elpat = re.compile(r"(Mo|W|S|Se|Hf|O)(?:[A-Z_]|\d|$)")
cand = [n for n in pots if re.search(r"MoS|WS|WSe|MoSe|Hf|_S_|_O_|_Mo|_W_", n)]
print("TMD_HFO2_CANDIDATES", len(cand))
for c in cand[:50]:
    print("   POT:", c)

print("--- sample of all potentials (first 60) ---")
for p in pots[:60]:
    print("   ", p)
print("PROBE_OK")
