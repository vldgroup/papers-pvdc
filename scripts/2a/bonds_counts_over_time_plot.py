import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

INPUT_CSV = "bonds_counts_over_time.csv"
OUTPUT_FIG = "bonds_counts_over_time.png"

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


df = pd.read_csv(INPUT_CSV)
df = df.sort_values("time_ps")

fig, ax = plt.subplots(figsize=(4.6, 2.75), dpi=200)

for i, stage in enumerate(STAGE_ORDER):
    stage_df = df[df["stage"] == stage]
    if stage_df.empty:
        continue

    t_axis = stage_df["time_ps"].to_numpy()

    if i == 0:
        cc_label = "C–C"
        ch_label = "C–H"
        ccl_label = "C–Cl"
    else:
        cc_label = None
        ch_label = None
        ccl_label = None

    ax.plot(t_axis, stage_df["cc_bonds"], color="#999999", label=cc_label, lw=1.5)
    ax.plot(t_axis, stage_df["ch_bonds"], color="#CC4668", label=ch_label, lw=2.0)
    ax.plot(t_axis, stage_df["ccl_bonds"], color="#34B091", label=ccl_label, lw=1.2)

ax.hlines(
    1764, 0, FULL_BOUNDS[-1], colors="#968C8C", linestyles="--", lw=1.0, alpha=0.8
)
for i in range(len(FULL_BOUNDS) - 1):
    if i % 2 == 1:
        ax.axvspan(
            FULL_BOUNDS[i], FULL_BOUNDS[i + 1], color="#968C8C", alpha=0.15, lw=0
        )

handles, labels = ax.get_legend_handles_labels()
order = [1, 0, 2]
ax.legend(
    [handles[idx] for idx in order],
    [labels[idx] for idx in order],
    loc="upper right",
    bbox_to_anchor=(0.98, 0.92),
    handletextpad=0.3,
    borderaxespad=0.2,
)

ax_ins = ax.inset_axes([0.5, 0.25, 0.3, 0.23])
for i, stage in enumerate(STAGE_ORDER):
    stage_df = df[df["stage"] == stage]
    if stage_df.empty:
        continue

    t_axis = stage_df["time_ps"].to_numpy()
    ax_ins.plot(t_axis, stage_df["cc_bonds"], color="#999999", lw=1.5)
    ax_ins.plot(t_axis, stage_df["ch_bonds"], color="#CC4668", lw=2.0)
    ax_ins.plot(t_axis, stage_df["ccl_bonds"], color="#34B091", lw=1.2)

for i in range(len(FULL_BOUNDS) - 1):
    if i % 2 == 1:
        ax_ins.axvspan(
            FULL_BOUNDS[i], FULL_BOUNDS[i + 1], color="#968C8C", alpha=0.15, lw=0
        )

ax_ins.set_xlim(30, 180)
ax_ins.set_ylim(1070, 1240)
ax_ins.yaxis.tick_right()
ax_ins.tick_params(labelsize=9.5, length=5)

ax.set_xlabel("Simulation time (ps)")
ax.set_ylabel("Number of bonds")
ax.set_xlim(0, FULL_BOUNDS[-1])
ax.set_ylim(-100, 1850)

mark_inset(ax, ax_ins, loc1=1, loc2=1, ec="#968C8C", lw=1.0, alpha=0.8, linestyle="--")
mark_inset(ax, ax_ins, loc1=3, loc2=3, ec="#968C8C", lw=1.0, alpha=0.8, linestyle="--")

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
