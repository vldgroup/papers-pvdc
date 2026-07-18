import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CSV = "reaction_type.csv"
OUTPUT_FIG = "reaction_type.png"

TIME_THRESHOLD = 178.76

plt.rcParams.update(
    {
        "axes.labelsize": 13,
        "axes.linewidth": 1.0,
        "axes.facecolor": "w",
        "axes.ymargin": 0.1,
        "axes.xmargin": 0.1,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 9.5,
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

count_data = np.zeros((2, 6), dtype=int)

for _, row in df.iterrows():
    time_ps = float(row["Time_ps"])
    connectivity = row["Connectivity"]
    euclidean_dist = float(row["C_C_euclidean_distance"])

    row_idx = 0 if time_ps < TIME_THRESHOLD else 1

    if connectivity == "1,2":
        col_idx = 1
    elif connectivity == "1,1":
        col_idx = 0
    elif connectivity == "1,3":
        col_idx = 2
    elif connectivity == "1,4":
        col_idx = 3
    elif connectivity == "others":
        if euclidean_dist < 6.7:
            col_idx = 4
        else:
            col_idx = 5
    else:
        continue

    count_data[row_idx, col_idx] += 1

fig, ax = plt.subplots(figsize=(3.9, 2.5), dpi=200)

plot_categories = [
    "1,2-elimination",
    "1,1/1,3/1,4-elimination",
    r"Others (proximal)",
    r"Others (distal)",
]
colors = ["#BD8BD2", "#04AC82", "#2D88D2", "#A6CEE3"]

grouped_widths = np.zeros((2, 4))
grouped_widths[:, 0] = count_data[:, 1]
grouped_widths[:, 1] = count_data[:, 0] + count_data[:, 2] + count_data[:, 3]
grouped_widths[:, 2] = count_data[:, 4]
grouped_widths[:, 3] = count_data[:, 5]

grouped_widths = grouped_widths[::-1]
x_labels = ["> 50%", "≤ 50%"]

totals = grouped_widths.sum(axis=1)
y_pos = np.array([0, 1])
left = np.zeros(2)

for i, cat in enumerate(plot_categories):
    widths = grouped_widths[:, i]
    bars = ax.barh(y_pos, widths, left=left, label=cat, color=colors[i], height=0.5)

    for j, rect in enumerate(bars):
        w = rect.get_width()
        if totals[j] == 0 or w == 0:
            continue
        percentage = (w / totals[j]) * 100

        if percentage > 8:
            ax.text(
                left[j] + w / 2,
                rect.get_y() + rect.get_height() / 2,
                f"{percentage:.1f}%",
                ha="center",
                va="center",
                color=("white" if i != 3 else "#333333"),
                fontsize=9.5,
                fontweight="bold",
            )
        elif percentage >= 1:
            ax.annotate(
                f"{percentage:.1f}%",
                xy=(left[j] + w / 2, rect.get_y() + rect.get_height() / 2 + 0.1),
                xytext=(
                    left[j] + w / 2 + 5,
                    rect.get_y() + rect.get_height() / 2 + 0.45,
                ),
                ha="left",
                va="center",
                fontsize=9.5,
                arrowprops=dict(
                    arrowstyle="-",
                    color="#333333",
                    lw=0.8,
                    connectionstyle="arc3,rad=0",
                ),
            )
    left += widths

ax.tick_params(axis="y", length=0)
ax.set_yticks(y_pos)
ax.set_yticklabels(x_labels)
ax.set_ylim(y_pos[0] - 0.5, y_pos[1] + 1.2)
ax.set_xlabel("Number of dehydrochlorination events")
ax.set_xlim(0, 760)

handles, labels = ax.get_legend_handles_labels()
order = [0, 2, 1, 3]
handles = [handles[i] for i in order]
labels = [labels[i] for i in order]

ax.legend(
    handles,
    labels,
    ncol=2,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.0),
    columnspacing=1,
    handletextpad=0.3,
    handlelength=1.0,
    labelspacing=0.2,
)

for i, total in enumerate(totals):
    ax.text(
        total + (max(totals) * 0.01),
        y_pos[i],
        f"{int(total)}",
        ha="left",
        va="center",
        color="#000000",
        fontsize=9.5,
    )

ax.set_yticks(y_pos)
ax.set_yticklabels(x_labels)

ax.text(
    -0.20,
    0.77,
    "C cluster\nformation",
    transform=ax.transAxes,
    ha="left",
    va="bottom",
    fontsize=9.5,
    multialignment="center",
)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
