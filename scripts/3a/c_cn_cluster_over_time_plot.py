import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CSV = "c_cn_cluster_over_time.csv"
OUTPUT_FIG = "c_cn_cluster_over_time.png"

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

FULL_BOUNDS = [0, 90.10, 271.90, 372.90, 645.60]
WINDOW_SIZE = 10


df = pd.read_csv(INPUT_CSV)
required_cols = {
    "stage",
    "time_ps",
    "largest_cluster_pct",
    "cn4_pct",
    "cn3_pct",
    "cn2_pct",
}
if not required_cols.issubset(df.columns):
    missing = sorted(required_cols - set(df.columns))
    raise ValueError(f"Input CSV is missing required columns: {missing}")

df = df.sort_values("time_ps")

fig, ax = plt.subplots(figsize=(4, 2.5), dpi=200)

for i, stage in enumerate(STAGE_ORDER):
    stage_df = df[df["stage"] == stage]
    if stage_df.empty:
        continue

    t_axis = stage_df["time_ps"].to_numpy()

    ax.plot(
        t_axis,
        stage_df["largest_cluster_pct"],
        color="#A0A0A0",
        lw=1.5,
        alpha=0.9,
        label="Largest cluster" if i == 0 else None,
    )

    for key, color, label in [
        ("cn4_pct", "#408AB3", r"CN$ = 4$"),
        ("cn3_pct", "#AD84C6", r"CN$ = 3$"),
        ("cn2_pct", "#A6B728", r"CN$ = 2$"),
    ]:
        data = stage_df[key].to_numpy()
        if len(data) >= WINDOW_SIZE:
            smoothed = np.convolve(
                data, np.ones(WINDOW_SIZE) / WINDOW_SIZE, mode="valid"
            )
            t_smooth = t_axis[: len(smoothed)]
        else:
            smoothed = data
            t_smooth = t_axis

        ax.plot(
            t_smooth,
            smoothed,
            color=color,
            lw=1.5,
            alpha=0.85,
            label=label if i == 0 else None,
        )

for i in range(len(FULL_BOUNDS) - 1):
    if i % 2 == 1:
        ax.axvspan(
            FULL_BOUNDS[i], FULL_BOUNDS[i + 1], color="#968C8C", alpha=0.12, lw=0
        )

ax.axhline(y=100, color="#968C8C", ls="--", lw=1.0, alpha=0.5)

handles, labels = ax.get_legend_handles_labels()
order = [1, 2, 3, 0]
ax.legend(
    [handles[idx] for idx in order],
    [labels[idx] for idx in order],
    bbox_to_anchor=(0.52, 0.74),
    handletextpad=0.3,
    borderaxespad=0.2,
    handlelength=1.8,
)

ax.set_xlabel("Simulation time (ps)")
ax.set_ylabel("Percentage C atoms (%)")
ax.set_xlim(0, FULL_BOUNDS[-1])
ax.set_ylim(-2, 105)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
