#!/usr/bin/env python3
"""Figure for the J0437 search: controls and injection recovery.

Both statistics are shown: the fixed-orientation filter and the F-statistic.
Recovery fractions are computed exactly from the saved sums.
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _common import HERE, ROOT, statistic_coefficients

RUN = ROOT / "results/j0437_search"
# Beside the article where the repository keeps articles in folders,
# otherwise in figures/.
ARTICLE_FIGURES = ROOT / "articles/j0437-one-orbital-template/figures"
OUT = (ARTICLE_FIGURES if ARTICLE_FIGURES.is_dir() else ROOT / "figures") \
    / "j0437_search_results.png"
INK, MUTED = "#0b0b0b", "#52514e"
SERIES = dict(fixed_orientation=(r"$2\mathcal{G}$", "#2a78d6"),
              f_statistic=(r"$2\mathcal{F}$", "#eb6834"))
CONTROLS = dict(fixed_orientation="fixed_orientation", f_statistic="two_f")

f_run = json.loads((RUN / "f_statistic.json").read_text())
summary = json.loads((HERE / "results/summary.json").read_text())
coefficients = statistic_coefficients(f_run)

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
                     "legend.fontsize": 8, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False})
figure, (left, right) = plt.subplots(1, 2, figsize=(7.25, 3.0),
                                     constrained_layout=True)

# (a) the 40 frequency controls and the target, one row per statistic
rng = np.random.default_rng(0)
for row, key in enumerate(SERIES):
    label, colour = SERIES[key]
    controls = np.array(f_run["controls"][CONTROLS[key]])
    target = summary[key]["target"]
    left.scatter(controls, row + rng.uniform(-0.16, 0.16, controls.size), s=16,
                 color=colour, alpha=0.75, linewidths=0,
                 label="40 frequency controls" if row == 0 else None)
    left.scatter([target], [row], s=62, marker="D", color=colour,
                 edgecolors=INK, linewidths=1.2, zorder=3,
                 label="target" if row == 0 else None)
    left.annotate(f"target {target:.2f}", (target, row), xytext=(0, 13),
                  textcoords="offset points", ha="center", color=INK, fontsize=8)
left.set_yticks([0, 1], [SERIES[k][0] for k in SERIES])
left.tick_params(axis="y", length=0, labelcolor=INK)
left.set_ylim(1.6, -0.6)
left.set_xlim(0, 18.5)
left.set_xlabel("statistic")
left.set_title("(a) Target and frequency controls", color=INK)
legend = left.legend(frameon=False, loc="lower right", handletextpad=0.3)
for handle in legend.legend_handles:
    handle.set_color(MUTED)

# (b) fraction of injections above the loudest control
grid = np.geomspace(1.5e-27, 1.6e-25, 500)
PLACE = dict(fixed_orientation=(4.5, 0.52, "left"), f_statistic=(0.7, 0.86, "right"))
for key, (label, colour) in SERIES.items():
    a, b, c = coefficients[key]
    threshold = summary[key]["loudest_control"]
    limits = summary[key]["h0_95"]
    fraction = np.array([np.mean(a * h * h + b * h + c > threshold) for h in grid])
    low, high = limits["above_loudest_control_bootstrap_90pct"]
    h95 = limits["above_loudest_control"]
    right.axvspan(low / 1e-26, high / 1e-26, color=colour, alpha=0.14, linewidth=0)
    right.plot(grid / 1e-26, fraction, color=colour, linewidth=2)
    right.plot([h95 / 1e-26] * 2, [0, 0.95], color=colour, linewidth=1,
               linestyle="--")
    x, y, align = PLACE[key]
    right.text(x, y, rf"{label}" "\n" rf"${h95 / 1e-26:.2f}\times10^{{-26}}$",
               color=INK, fontsize=8, ha=align, va="center", linespacing=1.3)
right.axhline(0.95, color=MUTED, linestyle=":", linewidth=1)
right.text(0.21, 0.965, "95%", color=MUTED, fontsize=8, va="bottom")
right.set_xscale("log")
right.set_xlim(0.2, 22)
right.set_ylim(0, 1.03)
right.set_xticks([0.2, 0.5, 1, 2, 5, 10], ["0.2", "0.5", "1", "2", "5", "10"])
right.minorticks_off()
right.set_xlabel(r"injected $h_0$ ($10^{-26}$)")
right.set_ylabel("fraction above loudest control")
right.set_title("(b) Recovery of 200 injections", color=INK)

OUT.parent.mkdir(parents=True, exist_ok=True)
figure.savefig(OUT, dpi=220, facecolor="white")
print(f"wrote {OUT.relative_to(ROOT)}")
