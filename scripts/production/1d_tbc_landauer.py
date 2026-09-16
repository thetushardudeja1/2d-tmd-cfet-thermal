"""PHASE 1d-A - proper Landauer thermal boundary conductance, WITH UNITS.

The earlier script printed 'LANDAUER_G_rel', which was NOT a physical TBC: it lacked
the hbar*omega weighting and the area normalisation. This fixes that.

FORMULA
    G = (1/2pi) * INT domega * hbar*omega * T(omega) * dn_BE/dT          [W/K]
With x = hbar*omega/(kB*T),  dn/dT = x e^x / (T (e^x-1)^2), and hbar*omega = x kB T:
    G = (kB/2pi) * INT domega * T(omega) * x^2 e^x / (e^x - 1)^2         [W/K]
Then TBC = G / A                                                         [W/m^2/K]

SANITY CHECK - the quantum of thermal conductance per ballistic mode:
    g0 = pi^2 kB^2 T / (3 h)  = 9.456e-13 * T  W/K   -> 2.837e-10 W/K at 300 K
A single perfectly transmitting mode must give exactly g0.

usage:  atkpython 1d_tbc_landauer.py <file.hdf5> [area_Ang2] [T]
"""
import sys
import numpy
from QuantumATK import *

HBAR = 1.054571817e-34      # J s
KB = 1.380649e-23           # J/K
QE = 1.602176634e-19        # C  (eV -> J)
H = 6.62607015e-34          # J s


def landauer_G(E_eV, T_of_E, temperature=300.0):
    """Landauer thermal conductance in W/K. E_eV in eV, T_of_E dimensionless."""
    E = numpy.asarray(E_eV, dtype=float) * QE          # J
    Tr = numpy.asarray(T_of_E, dtype=float)
    omega = E / HBAR                                   # rad/s
    x = E / (KB * temperature)
    # guard the x->0 limit: x^2 e^x/(e^x-1)^2 -> 1
    x = numpy.where(x < 1e-8, 1e-8, x)
    xc = numpy.clip(x, 0.0, 700.0)
    ex = numpy.exp(xc)
    w = (x ** 2) * ex / (ex - 1.0) ** 2                # dimensionless window
    integrand = KB * Tr * w                            # J/K
    return numpy.trapz(integrand, omega) / (2.0 * numpy.pi)   # W/K


def g0(temperature=300.0):
    """Quantum of thermal conductance per mode, W/K."""
    return (numpy.pi ** 2) * (KB ** 2) * temperature / (3.0 * H)


# ---------------- self-test: one perfectly transmitting acoustic-like mode -------
def self_test(temperature=300.0):
    # A single mode transmitting perfectly from 0 up to a high cutoff must give g0.
    E = numpy.linspace(0.0, 0.60, 20001)      # eV; well above kB*T (~0.026 eV)
    Tr = numpy.ones_like(E)
    G = landauer_G(E, Tr, temperature)
    ref = g0(temperature)
    print("SELFTEST  G(1 perfect mode) = %.4e W/K" % G)
    print("SELFTEST  g0 (quantum)      = %.4e W/K" % ref)
    print("SELFTEST  ratio             = %.4f   (must be ~1.000)" % (G / ref))
    return abs(G / ref - 1.0) < 0.02


if __name__ == "__main__":
    print("TBC_START")
    ok = self_test(300.0)
    print("SELFTEST_PASS" if ok else "SELFTEST_FAIL")

    if len(sys.argv) > 1:
        fname = sys.argv[1]
        area_A2 = float(sys.argv[2]) if len(sys.argv) > 2 else None
        temperature = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0

        pts = nlread(fname, PhononTransmissionSpectrum)[-1]
        Tr = numpy.array(pts.evaluate())
        E = numpy.array(pts.energies().inUnitsOf(eV))
        G = landauer_G(E, Tr, temperature)
        print("FILE %s" % fname)
        print("G_total   = %.6e W/K   (%.3f x g0)" % (G, G / g0(temperature)))
        if area_A2:
            A = area_A2 * 1e-20                      # Ang^2 -> m^2
            tbc = G / A
            print("AREA      = %.3f Ang^2 = %.4e m^2" % (area_A2, A))
            print("TBC       = %.4e W/m^2/K = %.2f MW/m^2/K" % (tbc, tbc / 1e6))
        else:
            print("(no area given -> TBC not computed; pass area in Ang^2 as arg 2)")
    print("TBC_DONE")
