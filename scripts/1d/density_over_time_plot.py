import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CSV = "density_over_time.csv"
OUTPUT_FIG = "density_over_time.png"

plt.rcParams.update(
    {
        "axes.labelsize": 13,
        "axes.linewidth": 1.0,
        "axes.facecolor": "w",
        "axes.ymargin": 0.1,
        "axes.xmargin": 0.1,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 10.5,
        "legend.fancybox": False,
        "legend.edgecolor": "k",
        "patch.linewidth": 0.8,
        "legend.borderaxespad": 0.8,
        "legend.title_fontsize": 12,
        "lines.linewidth": 0.8,
        "xtick.direction": "in",
        "xtick.major.size": 10,
        "xtick.major.width": 0.8,
        "xtick.minor.size": 8,
        "xtick.minor.width": 0.6,
        "ytick.direction": "in",
        "ytick.major.size": 10,
        "ytick.major.width": 0.8,
        "ytick.minor.size": 8,
        "ytick.minor.width": 0.6,
        "xtick.major.pad": 3,
        "ytick.major.pad": 3,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial"],
        "lines.markersize": 13,
    }
)

STAGE_ORDER = [
    "npt_ramp_1200",
    "nvt_ramp_3000",
    "nvt_hold_3000",
    "npt_quench_300",
]

WINDOW_SIZE = 5

df = pd.read_csv(INPUT_CSV)
df = df.sort_values("time_ps")

bounds = []
for stage in STAGE_ORDER:
    stage_df = df[df["stage"] == stage]
    if not stage_df.empty:
        bounds.append((stage_df["time_ps"].min(), stage_df["time_ps"].max()))

fig, ax = plt.subplots(figsize=(4, 2.5), dpi=200)

for stage in STAGE_ORDER:
    stage_df = df[df["stage"] == stage]
    if stage_df.empty:
        continue

    data = stage_df["density_g_cm3"].values
    time = stage_df["time_ps"].values

    if len(data) >= WINDOW_SIZE:
        smoothed_data = np.convolve(data, np.ones(WINDOW_SIZE), "valid") / WINDOW_SIZE
        time_smoothed = time[: len(smoothed_data)]
    else:
        smoothed_data = data
        time_smoothed = time

    ax.plot(
        time_smoothed,
        smoothed_data,
        color="#DF5372",
        lw=1.5,
        label="Simulation" if stage == STAGE_ORDER[0] else None,
    )

for i, (start, end) in enumerate(bounds):
    if i % 2 == 1:
        ax.axvspan(start, end, color="#968C8C", alpha=0.15, lw=0)

ax.set_xlabel("Simulation time (ps)")
ax.set_ylabel(r"Density (g cm$^{-3}$)")
ax.set_xlim(df["time_ps"].min(), df["time_ps"].max())
ax.set_ylim(1.13, 1.99)

ax.scatter(
    [40],
    [1.859],
    marker="h",
    s=40,
    color="#FFFFFF88",
    edgecolor="#000000",
    label="Exp. PVDC",
)

ax.scatter(
    [580] * 5,
    [1.19, 1.30, 1.32, 1.38, 1.59],
    marker="o",
    s=30,
    color="#FFFFFF88",
    edgecolor="#000000",
    label="Exp. hard C (1951)",
)

ax.scatter(
    [580] * 4,
    [1.43, 1.46, 1.59, 1.70],
    marker="^",
    s=35,
    color="#FFFFFF88",
    edgecolor="#000000",
    label="Exp. hard C (1964)",
)

ax.legend(
    loc="lower left",
    bbox_to_anchor=(0.04, 0.015),
    ncol=2,
    columnspacing=0.4,
    handletextpad=0.3,
    borderaxespad=0.2,
    handlelength=1.2,
)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
