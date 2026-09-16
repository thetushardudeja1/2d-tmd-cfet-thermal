"""Verify BOTH relaxed interfaces before running transmission.

nlsave APPENDS to HDF5, so mos2_hfo2_v10.hdf5 contains BOTH the v10 (Hf-terminated,
collapsed gap 0.81 A) and v11 (O-terminated, gap 4.70 A) structures. We must confirm the
file we carry forward is the O-terminated one, not the failed Hf-terminated one.

Checks per interface:
  * composition and atom count
  * TOP oxide termination (must be O)
  * vdW gap (oxide top -> lowest chalcogen) - must be ~4-5 A, NOT collapsed
  * TMD layer spread (must be ~2*dz, i.e. flat)
"""
import numpy
from QuantumATK import *

CASES = [
    ('iface_MoS2_HfO2.hdf5', 'MoS2', 'Mo', 'S',  3.129),
    ('WSe2_hfo2_v12.hdf5',   'WSe2', 'W',  'Se', 3.340),
]

print("VERIFY_START")
ok_all = True
for fn, tag, m_sym, x_sym, ideal in CASES:
    cs = nlread(fn, BulkConfiguration)
    print("=== %s : %d config(s) in file" % (fn, len(cs)))
    c = cs[-1]
    p = numpy.array(c.cartesianCoordinates().inUnitsOf(Angstrom))
    s = [str(e.symbol()) for e in c.elements()]
    comp = {e: s.count(e) for e in sorted(set(s))}

    ox_z = [p[i][2] for i, t in enumerate(s) if t in ("Hf", "O")]
    x_z = [p[i][2] for i, t in enumerate(s) if t == x_sym]
    tmd_z = [p[i][2] for i, t in enumerate(s) if t in (m_sym, x_sym)]
    z_ox = max(ox_z)
    gap = min(x_z) - z_ox
    spread = max(tmd_z) - min(tmd_z)

    # termination of the oxide top surface
    top = sorted(set(s[i] for i in range(len(s))
                     if p[i][2] > z_ox - 0.6 and s[i] in ("Hf", "O")))

    good = (top == ["O"]) and (3.0 < gap < 6.0) and (abs(spread - ideal) < 0.8)
    ok_all = ok_all and good
    print("    %-6s atoms=%d  %s" % (tag, len(s), comp))
    print("    oxide TOP termination : %s   %s" % ("/".join(top), "OK" if top == ["O"] else "WRONG"))
    print("    vdW gap               : %.3f A  %s" % (gap, "OK" if 3.0 < gap < 6.0 else "BAD"))
    print("    %-4s layer spread      : %.3f A (ideal %.3f)  %s"
          % (tag, spread, ideal, "OK" if abs(spread - ideal) < 0.8 else "DISTORTED"))
    print("    VERDICT: %s" % ("USABLE" if good else "NOT USABLE"))
    if good:
        out = 'FINAL_iface_%s_HfO2.hdf5' % tag
        nlsave(out, c)
        print("    -> saved clean copy: %s" % out)
    print()

print("ALL_INTERFACES_OK" if ok_all else "SOME_INTERFACES_BAD")
print("VERIFY_DONE")
