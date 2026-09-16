# Build publication-style figures from extracted QuantumATK data.
# Run locally: "C:/Users/TUSHAR/anaconda3/python.exe" scripts/production/make_figures.py
import csv
import os
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RAW = os.path.join(ROOT, "figures", "_raw_data")
OUT = os.path.join(ROOT, "figures")

MOS2_COLOR = "#2166ac"
WSE2_COLOR = "#b2182b"

plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def load_xy(path, xcol=0, ycol=1):
    with open(path) as f:
        r = list(csv.reader(f))
    rows = r[1:]
    data = np.array(rows, dtype=float)
    return data[:, xcol], data[:, ycol]


def fig_dos():
    eM, dM = load_xy(os.path.join(RAW, "MoS2_dos.csv"))
    eW, dW = load_xy(os.path.join(RAW, "WSe2_dos.csv"))
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(eM, dM, color=MOS2_COLOR, lw=1.6, label="MoS$_2$")
    ax.plot(eW, dW, color=WSE2_COLOR, lw=1.6, ls="--", label="WSe$_2$")
    ax.set_xlabel("Phonon frequency (cm$^{-1}$)")
    ax.set_ylabel("Phonon DOS (arb. units)")
    ax.set_title("Phonon density of states: MoS$_2$ vs WSe$_2$")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_phonon_dos.png"))
    plt.close(fig)
    print("fig1_phonon_dos.png done")


def fig_transmission():
    eM, tM = load_xy(os.path.join(RAW, "MoS2_HfO2_transmission_NTMD4.csv"))
    eW, tW = load_xy(os.path.join(RAW, "WSe2_HfO2_transmission_NTMD4.csv"))
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(eM, tM, color=MOS2_COLOR, lw=1.6, label="MoS$_2$/HfO$_2$ (TBC = 14.95 MW m$^{-2}$K$^{-1}$)")
    ax.plot(eW, tW, color=WSE2_COLOR, lw=1.6, ls="--", label="WSe$_2$/HfO$_2$ (TBC = 12.85 MW m$^{-2}$K$^{-1}$)")
    ax.set_xlabel("Phonon energy (cm$^{-1}$)")
    ax.set_ylabel("Transmission $T(E)$ (q-averaged)")
    ax.set_title("Ballistic phonon transmission spectrum at matched thickness\n(N$_{TMD}$ = 4 free layers, 6.3 nm gate stack)")
    ax.set_xlim(0, max(eM.max(), eW.max()))
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_transmission_spectrum.png"))
    plt.close(fig)
    print("fig2_transmission_spectrum.png done")


def fig_tbc_ratio():
    n_tmd = np.array([2, 4])
    mos2 = np.array([8.20, 14.95])
    wse2 = np.array([7.24, 12.85])
    ratio = mos2 / wse2

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.2))

    width = 0.35
    x = np.arange(len(n_tmd))
    ax1.bar(x - width / 2, mos2, width, color=MOS2_COLOR, label="MoS$_2$/HfO$_2$")
    ax1.bar(x + width / 2, wse2, width, color=WSE2_COLOR, label="WSe$_2$/HfO$_2$")
    ax1.set_xticks(x)
    ax1.set_xticklabels(["N$_{TMD}$=2\n(4 layers)", "N$_{TMD}$=4\n(6 layers)"])
    ax1.set_ylabel("TBC (MW m$^{-2}$K$^{-1}$)")
    ax1.set_title("Absolute TBC\n(still rising with TMD thickness)")
    ax1.legend(frameon=False, fontsize=9)

    ax2.plot(n_tmd, ratio, "o-", color="0.2", lw=2, ms=8)
    for xi, yi in zip(n_tmd, ratio):
        ax2.annotate("%.3f" % yi, (xi, yi), textcoords="offset points",
                     xytext=(0, 10), ha="center", fontsize=10)
    ax2.set_xticks(n_tmd)
    ax2.set_xlabel("N$_{TMD}$ (free TMD layers in cell)")
    ax2.set_ylabel("Ratio MoS$_2$ : WSe$_2$ TBC")
    ax2.set_title("Ratio is stable to 2.7%\nwhile absolutes move ~80%")
    ax2.set_ylim(1.0, 1.3)

    fig.suptitle("Material-driven TBC ratio: the paper-ready result", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_tbc_ratio.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig3_tbc_ratio.png done")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig_dos()
    fig_transmission()
    fig_tbc_ratio()
    print("ALL_FIGURES_DONE")
