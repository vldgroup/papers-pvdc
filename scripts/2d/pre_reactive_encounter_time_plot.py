from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CSV = "pre_reactive_encounter_time.csv"
OUTPUT_FIG = "pre_reactive_encounter_time.png"

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

np.random.seed(0)

df = pd.read_csv(INPUT_CSV)

df["1st_iso_time_ps"] = pd.to_numeric(df["1st_iso_time_ps"], errors="coerce")
df = df.dropna(subset=["Symbol", "1st_iso_time_ps"])
df = df[df["Symbol"].isin(["H", "Cl"])].reset_index(drop=True)

df_long = df.rename(columns={"1st_iso_time_ps": "Time_ps"}).copy()
H_1st = df_long[df_long["Symbol"] == "H"]["Time_ps"].to_numpy()
Cl_1st = df_long[df_long["Symbol"] == "Cl"]["Time_ps"].to_numpy()

H_1st = H_1st[np.isfinite(H_1st) & (H_1st > 0)]
Cl_1st = Cl_1st[np.isfinite(Cl_1st) & (Cl_1st > 0)]

fig, ax = plt.subplots(figsize=(3.7, 2.5), dpi=200)

category_centers = [1, 2]
scatter_positions = [0.82, 1.82]
box_positions = [1.18, 2.18]

data = [H_1st, Cl_1st]

box = ax.boxplot(
    data,
    positions=box_positions,
    widths=0.3,
    patch_artist=True,
    showfliers=False,
    showmeans=True,
    meanline=True,
    meanprops=dict(color="black", linestyle=":", linewidth=1.1),
    medianprops=dict(color="black", linewidth=1.1),
    boxprops=dict(linewidth=1),
    whiskerprops=dict(linewidth=1),
    capprops=dict(linewidth=1),
)

colors = ["#E5B5B8", "#8DD3C7"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.85)


def jitter_scatter(values, xpos, color):
    if len(values) == 0:
        return
    x = np.random.normal(xpos, 0.04, size=len(values))
    ax.scatter(x, values, s=8, color=color, alpha=0.35, edgecolor="none")


jitter_scatter(H_1st, scatter_positions[0], "#C1666B")
jitter_scatter(Cl_1st, scatter_positions[1], "#264653")

ax.set_yscale("log")
ax.set_xlim(0.3, 2.7)
ax.set_ylim(0.006, 15)
ax.set_xticks(category_centers)
ax.set_xticklabels(["H", "Cl"])
ax.set_xlabel("Radical species")
ax.set_ylabel("PRE time (ps)")
ax.grid(axis="y", linestyle="--", alpha=0.7)
ax.xaxis.grid(False)

for spine in ax.spines.values():
    spine.set_color("black")

line_legend_elements = [
    Line2D([0], [0], color="black", linestyle="-", linewidth=1, label="Median"),
    Line2D([0], [0], color="black", linestyle=":", linewidth=1, label="Mean"),
]
ax.legend(handles=line_legend_elements, loc="upper left", bbox_to_anchor=(0.05, 0.97), handletextpad=0.3, borderaxespad=0.2)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
