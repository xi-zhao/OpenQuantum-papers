"""Publication figures generated only from the verified analysis files."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spectral_analysis import constants

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT / "figs"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
MARKERS = ["o", "s", "^"]
plt.rcParams.update({
  "font.family": "DejaVu Sans", "font.size": 8,
  "axes.labelsize": 8, "axes.titlesize": 8.5,
  "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
  "axes.spines.top": False, "axes.spines.right": False,
  "axes.linewidth": .6, "axes.edgecolor": "#888888",
  "text.color": "#222222", "axes.labelcolor": "#222222",
  "xtick.color": "#52514e", "ytick.color": "#52514e",
  "lines.linewidth": 1.2, "lines.markersize": 4.5,
  "pdf.fonttype": 42, "ps.fonttype": 42,
  "savefig.facecolor": "white",
})


def save(fig, name):
  fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight")
  fig.savefig(FIGS / f"{name}.png", dpi=200, bbox_inches="tight")
  plt.close(fig)


def splitting():
  fig, ax = plt.subplots(figsize=(3.3, 2.4), layout="constrained")
  yroots = np.roots([1, -504, 27504])
  angles = 2 * np.pi * np.arange(3) / 3
  exact = np.concatenate([r ** (1/3) * np.exp(1j * angles) for r in yroots])
  truncated = 504 ** (1/3) * np.exp(1j * angles)
  ax.scatter(exact.real, exact.imag, c=COLORS[0], s=28, marker="o",
             edgecolors="white", linewidths=.6, label="Exact: two triplets", zorder=3)
  ax.scatter(truncated.real, truncated.imag, c=COLORS[1], s=32, marker="x",
             linewidths=1.2, label="Linear-in-z truncation", zorder=4)
  ax.scatter([0], [0], c="#444444", s=24, marker="D", edgecolors="white", linewidths=.5)
  ax.annotate("zero: multiplicity 2\n(truncation: 5)", (0, 0), xytext=(.6, -2.9),
              fontsize=7, arrowprops={"arrowstyle": "-", "color": "#999999", "lw": .6})
  ax.axhline(0, color="#dddddd", lw=.5, zorder=0)
  ax.axvline(0, color="#dddddd", lw=.5, zorder=0)
  ax.set(xlabel=r"Re $[\mu/z^{1/3}]$", ylabel=r"Im $[\mu/z^{1/3}]$",
         title=r"$n=7,\ d=2$: all retained orders matter", ylim=(-8, 9.3), xlim=(-5.5, 10.5))
  ax.legend(frameon=False, loc="upper left", handletextpad=.4)
  ax.set_aspect("equal", adjustable="box")
  save(fig, "splitting")
  rows = [("exact", a.real, a.imag) for a in exact]
  rows += [("truncation", a.real, a.imag) for a in truncated]
  with (ROOT / "data" / "splitting_table.csv").open("w") as f:
    writer = csv.writer(f, lineterminator="\n")
    writer.writerow(["spectrum", "real_scaled_mu", "imag_scaled_mu"])
    writer.writerows(rows)
    writer.writerow(["exact_zero_multiplicity_2", 0, 0])
    writer.writerow(["truncated_zero_multiplicity_5", 0, 0])


def scaling():
  perron = json.loads((ROOT / "data" / "verified_perron.json").read_text())
  precision = json.loads((ROOT / "data" / "verified_precision.json").read_text())
  fig, axes = plt.subplots(1, 3, figsize=(7, 2.38), layout="constrained")
  table = []
  for idx, d in enumerate([2, 3, 4]):
    c = constants(d)
    rows = [r for r in perron["points"] if r["d"] == d]
    n = np.array([r["n"] for r in rows])
    mu = np.array([r["value"] for r in rows])
    root_ratio = mu / (c["C"] * n ** c["alpha"])
    scaled_correction = n * (root_ratio - 1)
    # Distinct markers and a companion CSV preserve information in grayscale.
    style = dict(color=COLORS[idx], marker=MARKERS[idx], markeredgecolor="white",
                 markeredgewidth=.5, label=rf"$d={d}$")
    axes[0].plot(1/n, root_ratio, **style)
    axes[1].plot(1/n, scaled_correction, **style)
    correction_limit = c["first_correction"] / c["C"]
    axes[1].plot([0, .064], [correction_limit]*2, color=COLORS[idx], lw=.8, ls="--")
    axes[1].text(.049, scaled_correction[0]+.025, rf"$d={d}$", fontsize=7)
    pr = [r for r in precision["points"] if r["d"] == d]
    pn = np.array([r["n"] for r in pr])
    gap = np.array([r["gap"] for r in pr])
    axes[2].plot(1/pn, gap / (c["internal_gap"] * pn ** c["beta"]), **style)
    for r, rr, cor in zip(rows, root_ratio, scaled_correction):
      table.append([d, r["n"], r["value"], r["lower"], r["upper"], rr, cor,
                    correction_limit, r["gap_one"], r["gap_two"]])
  axes[0].axhline(1, color="#888888", ls="--", lw=.8)
  axes[2].axhline(1, color="#888888", ls="--", lw=.8)
  axes[0].set(title="(a) Leading spectral edge", ylabel=r"$\kappa_{n,d}/(C_d n^{\alpha_d})$")
  axes[1].set(title="(b) First finite-size correction", ylabel=r"$n[\kappa_{n,d}/(C_d n^{\alpha_d})-1]$")
  axes[2].set(title="(c) Symmetric-block spacing", ylabel=r"$\delta_{n,d}/(2a_d\Omega_d n^{\beta_d})$")
  for ax in axes:
    ax.set_xlabel(r"$1/n$")
    ax.set_xlim(-.001, .066)
    ax.set_xticks([0, .02, .04, .06])
    ax.tick_params(length=3, width=.6)
  axes[0].legend(frameon=False, loc="lower left")
  axes[2].legend(frameon=False, loc="upper left")
  axes[0].text(.056, .984, r"$d=2$", fontsize=7)
  axes[0].text(.056, .966, r"$d=3$", fontsize=7)
  axes[0].text(.056, .946, r"$d=4$", fontsize=7)
  # Panel (c) curves converge; legend supplies identity without detached labels.
  save(fig, "scaling")
  with (ROOT / "data" / "scaling_table.csv").open("w") as f:
    writer = csv.writer(f, lineterminator="\n")
    writer.writerow(["d", "n", "perron", "lower", "upper", "leading_ratio",
                     "scaled_correction", "predicted_correction", "gap_one", "gap_two"])
    writer.writerows(table)


if __name__ == "__main__":
  FIGS.mkdir(exist_ok=True)
  splitting()
  scaling()
