import matplotlib.pyplot as plt
import pandas as pd

INPUT_CSV = "md_scheme_with_temp.csv"
OUTPUT_FIG = "md_scheme_with_temp.png"

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

df = pd.read_csv(INPUT_CSV)
observed = df[df["stage"] != "target"].sort_values("time_ps")
target = df[df["stage"] == "target"].sort_values("time_ps")

fig, ax = plt.subplots(figsize=(4, 2.5), dpi=200)

ax.plot(
    observed["time_ps"],
    observed["temp_k"],
    color="#999999",
    linewidth=1.2,
    alpha=0.75,
    label="Observed temperature",
)
ax.plot(
    target["time_ps"],
    target["temp_k"],
    color="#DF5372",
    linewidth=1.5,
    label="Target temperature",
)

ax.text(45, 3350, "NPT", ha="center", va="bottom", fontsize=11)
ax.axvspan(90, 270, color="#968C8C", alpha=0.15, lw=0)
ax.text(180, 3350, "NVT*", ha="center", va="bottom", fontsize=11)
ax.text(320, 3350, "NVT*", ha="center", va="bottom", fontsize=11)
ax.axvspan(370, 640, color="#968C8C", alpha=0.15, lw=0)
ax.text(505, 3350, "NPT*", ha="center", va="bottom", fontsize=11)

ax.set_xlabel("Simulation time (ps)")
ax.set_ylabel("Temperature (K)")
ax.set_xlim(0, 640)
ax.set_ylim(0, 3800)
ax.set_xticks([90, 270, 370, 640])
ax.set_yticks([0, 300, 1200, 3000])

ax.grid(True, linestyle="--", alpha=0.8, axis="y", lw=1.0, color="#968C8C")
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[::-1],
    labels[::-1],
    loc="lower center",
    bbox_to_anchor=(0.5, 0.03),
    labelspacing=0.3,
    borderpad=0.3,
)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
