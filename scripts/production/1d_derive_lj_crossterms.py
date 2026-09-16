"""Derive LJ cross-terms for MoS2/HfO2 and WSe2/HfO2 interfaces (Phase 1d).

METHOD (exactly D1 = Suryavanshi/Pop, JAP 126 055107 (2019), supplementary S2):
    Universal Force Field (UFF, Rappe et al., JACS 114, 10024 (1992))
    + Lorentz-Berthelot mixing:
        eps(A-B) = sqrt( eps(A) * eps(B) )        [geometric mean]
        sigma(A-B) = ( sigma(A) + sigma(B) ) / 2  [arithmetic mean]

UFF tabulates x_I (the vdW MINIMUM position, r_min) and D_I (well depth, kcal/mol).
The 12-6 LJ sigma is therefore:
        sigma = x / 2**(1/6)

VALIDATION: this script first reproduces D1's published Table S1. If our UFF inputs
are right, Mo-Mo, S-S, Mo-O, S-O, Mo-Si, S-Si must come out matching D1. Only then are
the NEW pairs (Mo-Hf, S-Hf, W-*, Se-*) trustworthy.

Run locally: python 1d_derive_lj_crossterms.py
"""
import math

KCAL_PER_MOL_TO_MEV = 43.3641          # 1 kcal/mol in meV
TWO_POW_1_6 = 2.0 ** (1.0 / 6.0)       # 1.122462

# --- UFF nonbond parameters: element -> (x_I in Angstrom, D_I in kcal/mol) ----------
# NOTE: these are read from the UFF table (Rappe 1992). They are VALIDATED below
# against D1's published cross-terms; if a value were wrong the validation would fail.
UFF = {
    "Mo": (3.052, 0.056),   # Mo6+6
    "S":  (4.035, 0.274),   # S_3+2
    "Hf": (2.798, 0.072),   # Hf3+4
    "O":  (3.500, 0.060),   # O_3
    "Si": (4.295, 0.402),   # Si3
    "W":  (2.734, 0.067),   # W_6+6
    "Se": (4.205, 0.291),   # Se3+2
    "Al": (4.499, 0.505),   # Al3
    "N":  (3.660, 0.069),   # N_3
}


def sigma_eps(el):
    """UFF (x, D) -> LJ (sigma in Ang, epsilon in meV)."""
    x, D = UFF[el]
    return x / TWO_POW_1_6, D * KCAL_PER_MOL_TO_MEV


def mix(a, b):
    """Lorentz-Berthelot mixing of two elements."""
    sa, ea = sigma_eps(a)
    sb, eb = sigma_eps(b)
    return (sa + sb) / 2.0, math.sqrt(ea * eb)


# --- 1. VALIDATION against D1 Table S1 --------------------------------------------
# (pair, D1 epsilon meV, D1 sigma Ang)
D1_TABLE = [
    ("Mo", "Si", 6.52, 3.274),
    ("Mo", "O",  2.51, 2.920),
    ("S",  "Si", 14.39, 3.721),
    ("S",  "O",  5.58, 3.385),
    ("Mo", "Al", 7.29, 3.362),
    ("Mo", "N",  2.69, 2.991),
    ("S",  "Al", 16.13, 3.803),
    ("S",  "N",  5.96, 3.428),
    ("Mo", "Mo", 2.40, 2.720),
    ("S",  "S",  11.88, 5.396),   # <-- sigma looks like a TYPO in D1 (see note below)
    ("Mo", "S",  5.34, 3.544),    # <-- also inconsistent; unused (intra-MoS2 = SW)
]

print("=" * 78)
print("VALIDATION - reproduce D1 (Suryavanshi/Pop 2019) Table S1")
print("=" * 78)
print("%-8s %10s %10s %8s   %10s %10s %8s" %
      ("pair", "eps_ours", "eps_D1", "d%", "sig_ours", "sig_D1", "d%"))
for a, b, e_d1, s_d1 in D1_TABLE:
    s, e = mix(a, b)
    de = 100.0 * (e - e_d1) / e_d1
    ds = 100.0 * (s - s_d1) / s_d1
    flag = "" if (abs(de) < 3 and abs(ds) < 3) else "   <-- MISMATCH"
    print("%-8s %10.3f %10.3f %+7.1f   %10.3f %10.3f %+7.1f%s"
          % (a + "-" + b, e, e_d1, de, s, s_d1, ds, flag))

print("""
NOTES ON THE TWO MISMATCHES
  S-S   : our sigma = 3.595 A vs D1's printed 5.396 A. Our value is consistent with
          ALL FOUR of D1's own S-X cross-terms (S-O, S-Si, S-Al, S-N), each of which
          back-solves to sigma_S = 3.60 +/- 0.03 A. D1's '5.396' is therefore almost
          certainly a TYPO for 3.596. Our epsilon matches D1 exactly (11.88 meV).
  Mo-S  : D1's printed sigma (3.544) is not consistent with their own Mo and S values
          either. IRRELEVANT for us: Mo-S is an INTRA-MoS2 pair, handled by the
          Stillinger-Weber potential, never by an LJ cross-term.
  Every pair that we actually USE as a cross-term reproduces D1 to <3%.
""")

# --- 2. THE NEW CROSS-TERMS WE NEED -----------------------------------------------
print("=" * 78)
print("DERIVED CROSS-TERMS  (12-6 LJ:  V(r) = 4*chi*eps*[(sig/r)^12 - (sig/r)^6])")
print("=" * 78)

print("\n--- MoS2 / HfO2 interface ---")
for a, b in [("Mo", "Hf"), ("Mo", "O"), ("S", "Hf"), ("S", "O")]:
    s, e = mix(a, b)
    src = "D1-published" if b == "O" else "NEW (this work)"
    print("  %-6s eps = %8.4f meV   sigma = %7.4f Ang    [%s]" % (a + "-" + b, e, s, src))

print("\n--- WSe2 / HfO2 interface (for the p-tier) ---")
for a, b in [("W", "Hf"), ("W", "O"), ("Se", "Hf"), ("Se", "O")]:
    s, e = mix(a, b)
    print("  %-6s eps = %8.4f meV   sigma = %7.4f Ang    [NEW (this work)]" % (a + "-" + b, e, s))

print("\n--- MoS2 / WSe2 interface (tier-to-tier vdW gap) ---")
for a, b in [("Mo", "W"), ("Mo", "Se"), ("S", "W"), ("S", "Se")]:
    s, e = mix(a, b)
    print("  %-6s eps = %8.4f meV   sigma = %7.4f Ang    [NEW (this work)]" % (a + "-" + b, e, s))

print("\n--- element self-parameters (for reference) ---")
for el in ["Mo", "S", "W", "Se", "Hf", "O"]:
    s, e = sigma_eps(el)
    print("  %-3s sigma = %7.4f Ang   eps = %8.4f meV" % (el, s, e))

print("""
REMINDER: chi (interaction-strength scaling) MUST be swept 1 -> 4 as D1 does.
UFF mixing rules ignore polarization and are known to UNDER-estimate the real
coupling: D1 found experiment sits between chi=1 and chi=2.
""")
