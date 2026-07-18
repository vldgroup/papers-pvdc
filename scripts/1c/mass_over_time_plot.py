import matplotlib.pyplot as plt
import pandas as pd

INPUT_CSV = "mass_over_time.csv"
OUTPUT_FIG = "mass_over_time.png"

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

C_ATOM_COUNT = 3528 // 3
C_MASS_DA = 12.01

df = pd.read_csv(INPUT_CSV)
df["total_mass_kda"] = df["total_mass_da"] / 1000.0
df = df.sort_values("time_ps")

bounds = []
for stage in STAGE_ORDER:
    stage_df = df[df["stage"] == stage]
    if not stage_df.empty:
        bounds.append((stage_df["time_ps"].min(), stage_df["time_ps"].max()))

plt.figure(figsize=(4, 2.5), dpi=200)
plt.plot(
    df["time_ps"], df["total_mass_kda"], color="#DF5372", lw=1.5, label="Total mass"
)

for i, (start, end) in enumerate(bounds):
    if i % 2 == 1:
        plt.axvspan(start, end, color="#968C8C", alpha=0.15, lw=0)

c_only_mass_kda = (C_ATOM_COUNT * C_MASS_DA) / 1000.0
plt.hlines(
    c_only_mass_kda,
    xmin=df["time_ps"].min(),
    xmax=df["time_ps"].max(),
    colors="#968C8C",
    linestyles="--",
    lw=1.0,
    alpha=0.8,
    label="Complete reaction",
)

plt.xlabel("Simulation time (ps)")
plt.ylabel("Mass (kDa)")
plt.xlim(df["time_ps"].min(), df["time_ps"].max())
plt.ylim(0, 65.0)
plt.legend(
    loc="lower left",
    ncol=2,
    bbox_to_anchor=(0.08, 0.045),
    borderaxespad=0.1,
    columnspacing=0.8,
    handletextpad=0.4,
    handlelength=1.5,
)
plt.tight_layout()

plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
