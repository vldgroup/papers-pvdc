import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CSV = "event_based_lifetime_freq.csv"
OUTPUT_FIG = "event_based_lifetime_freq.png"

WINDOW = 200

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

events = pd.read_csv(INPUT_CSV)

events["Duration_ps"] = pd.to_numeric(events["Duration_ps"], errors="coerce")
events["Start_time_ps"] = pd.to_numeric(events["Start_time_ps"], errors="coerce")
events = events.dropna(subset=["Symbol", "Duration_ps", "Start_time_ps"])

# apply same persistence threshold as original notebook
PERSISTENCE_THRESHOLD_PS = 0.01
events = events[events["Duration_ps"] > PERSISTENCE_THRESHOLD_PS]

df_cl = events[events["Symbol"].astype(str).str.strip() == "Cl"].copy()
df_cl = df_cl.rename(columns={"Duration_ps": "lifetime", "Start_time_ps": "birth_time"})
df_cl = df_cl.sort_values("birth_time").reset_index(drop=True)

if len(df_cl) >= WINDOW:
    df_cl["rolling_mean"] = df_cl["lifetime"].rolling(window=WINDOW).mean()
else:
    df_cl["rolling_mean"] = df_cl["lifetime"].expanding().mean()

df_cl["event_freq"] = np.nan
if len(df_cl) > WINDOW:
    dt = df_cl["birth_time"].diff(WINDOW)
    df_cl["event_freq"] = WINDOW / dt.replace(0, np.nan)

fig, ax1 = plt.subplots(figsize=(4.3, 2.5), dpi=200)

mean_mask = df_cl["birth_time"] >= 95.0
df_mean_plot = df_cl.loc[mean_mask].copy()

if not df_mean_plot.empty:
    ax1.plot(
        df_mean_plot["birth_time"],
        df_mean_plot["rolling_mean"],
        color="#BD8BD2",
        lw=1.2,
        label="Mean lifetime",
    )

ax1.set_yscale("log")
ax1.set_xlim(90, 370)
ax1.set_ylim(0.01, 20)
ax1.set_xlabel("Simulation time (ps)")
ax1.set_ylabel("Mean Cl· lifetime (ps)")

ax2 = ax1.twinx()
ax2.plot(
    df_cl["birth_time"],
    df_cl["event_freq"],
    color="#1A7A82",
    lw=1.2,
    label="Formation rate",
)

ax2.set_ylabel("Cl· formation rate (ps$^{-1}$)")
ax2.set_ylim(-20, 320)

ax2.axvspan(149.84, 185.17, color="#968C8C", alpha=0.25, lw=0)
ax2.axvline(
    x=178.76,
    color="#968C8C",
    linestyle="--",
    lw=1.1,
    alpha=0.8,
    label="50% C clustering",
)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax1.legend(
    lines1 + lines2,
    labels1 + labels2,
    loc="upper right",
    bbox_to_anchor=(0.97, 0.99),
    handletextpad=0.3,
    borderaxespad=0.2,
    handlelength=1.5,
)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
print(f"Saved figure to {OUTPUT_FIG}")
