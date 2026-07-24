#!/usr/bin/env python3

"""
scatterplot.py

Reads ~/Stenophage/Summary/phage_summary.csv and makes:
1. Genome size vs GC%, colored by tRNA count
2. Genome size vs CDS count

Adds best-fit line, equation, R², and Pearson r.
Outputs plots to ~/Stenophage/Summary.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

summary_folder = Path("~/Stenophage/Summary").expanduser()
csv_file = summary_folder / "phage_summary.csv"

df = pd.read_csv(csv_file)

for col in ["genome_length_bp", "GC_percent", "tRNA_count", "CDS_count"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

def thousands(x, pos):
    return f"{int(x):,}"

def add_fit(ax, x, y, labels=None, eq_x=0.04, eq_y=0.95, ha="left", multialignment=None):
    m, b = np.polyfit(x, y, 1)
    y_pred = m * x + b
    r = np.corrcoef(x, y)[0, 1]
    r2 = r ** 2
    xs = np.linspace(x.min(), x.max(), 200)
    ax.plot(xs, m * xs + b, linewidth=2, color="black", linestyle="--")
    eq = f"y = {m:.3e}x + {b:.2f}\nR² = {r2:.3f}\nr = {r:.3f}"
    ax.text(eq_x, eq_y, eq, transform=ax.transAxes, va="top", ha=ha, fontsize=10, multialignment=multialignment, bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.85))

    if labels is not None:
        residuals = y - y_pred
        std = residuals.std()
        if std > 0:
            outlier_mask = residuals.abs() > 2 * std
            for xi, yi, name in zip(x[outlier_mask], y[outlier_mask], labels[outlier_mask]):
                ax.annotate(
                    str(name),
                    xy=(xi, yi),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=7,
                    color="black",
                )

df1 = df.dropna(subset=["genome_length_bp", "GC_percent", "tRNA_count", "file"])
df2 = df.dropna(subset=["genome_length_bp", "CDS_count", "file"])

fig, ax = plt.subplots(figsize=(7, 5.5))
sc = ax.scatter(df1["genome_length_bp"], df1["GC_percent"], c=df1["tRNA_count"], s=65, alpha=0.9, edgecolors="black", linewidths=0.5, cmap="viridis", vmin=0, vmax=15)
add_fit(ax, df1["genome_length_bp"], df1["GC_percent"], labels=df1["file"], eq_x=0.97, eq_y=0.97, ha="right", multialignment="left")
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("tRNA count")
cbar.set_ticks([0, 2, 4, 6, 8, 10, 12, 14])
ax.set_xlabel("Genome size (bp)")
ax.set_ylabel("GC content (%)")
ax.set_title("A Genome size vs GC content")
ax.xaxis.set_major_formatter(FuncFormatter(thousands))
fig.tight_layout()
fig.savefig(summary_folder / "scatter_genome_size_vs_gc_trna.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.scatter(df2["genome_length_bp"], df2["CDS_count"], s=65, alpha=0.9, edgecolors="black", linewidths=0.5)
add_fit(ax, df2["genome_length_bp"], df2["CDS_count"], labels=df2["file"], eq_x=0.04, eq_y=0.97, ha="left")
ax.set_xlabel("Genome size (bp)")
ax.set_ylabel("CDS count")
ax.set_title("B  Genome size vs CDS count")
ax.xaxis.set_major_formatter(FuncFormatter(thousands))
fig.tight_layout()
fig.savefig(summary_folder / "scatter_genome_size_vs_cds.png", dpi=300)
plt.close(fig)

df3 = df.dropna(subset=["genome_length_bp", "tRNA_count", "file"])

def plot_c(ax):
    ax.scatter(df3["genome_length_bp"], df3["tRNA_count"], s=65, alpha=0.9, edgecolors="black", linewidths=0.5, color="tab:blue")
    add_fit(ax, df3["genome_length_bp"], df3["tRNA_count"], labels=df3["file"], eq_x=0.04, eq_y=0.97, ha="left")
    ax.set_xlabel("Genome size (bp)")
    ax.set_ylabel("tRNA count")
    ax.set_title("C  Genome size vs tRNA count")
    ax.xaxis.set_major_formatter(FuncFormatter(thousands))
    ax.set_ylim(bottom=-1)
    yticks = [t for t in ax.get_yticks() if t >= 0]
    ax.set_yticks(yticks)

fig, ax = plt.subplots(figsize=(7, 5.5))
plot_c(ax)
fig.tight_layout()
fig.savefig(summary_folder / "scatter_genome_size_vs_trna.png", dpi=300)
plt.close(fig)

# -------------------------------------------------------------------
# Combined figure: two graphs on top, one on the bottom
# -------------------------------------------------------------------
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(2, 4)

ax = fig.add_subplot(gs[0, 0:2])
sc = ax.scatter(df1["genome_length_bp"], df1["GC_percent"], c=df1["tRNA_count"], s=65, alpha=0.9, edgecolors="black", linewidths=0.5, cmap="viridis", vmin=0, vmax=15)
add_fit(ax, df1["genome_length_bp"], df1["GC_percent"], labels=df1["file"], eq_x=0.97, eq_y=0.97, ha="right", multialignment="left")
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("tRNA count")
cbar.set_ticks([0, 2, 4, 6, 8, 10, 12, 14])
ax.set_xlabel("Genome size (bp)")
ax.set_ylabel("GC content (%)")
ax.set_title("A Genome size vs GC content")
ax.xaxis.set_major_formatter(FuncFormatter(thousands))

ax = fig.add_subplot(gs[0, 2:4])
ax.scatter(df2["genome_length_bp"], df2["CDS_count"], s=65, alpha=0.9, edgecolors="black", linewidths=0.5)
add_fit(ax, df2["genome_length_bp"], df2["CDS_count"], labels=df2["file"], eq_x=0.04, eq_y=0.97, ha="left")
ax.set_xlabel("Genome size (bp)")
ax.set_ylabel("CDS count")
ax.set_title("B  Genome size vs CDS count")
ax.xaxis.set_major_formatter(FuncFormatter(thousands))

ax = fig.add_subplot(gs[1, 1:3])
plot_c(ax)

fig.tight_layout()
fig.savefig(summary_folder / "scatter_genome_size_combined.png", dpi=300)
plt.close(fig)

print("Saved:")
print(summary_folder / "scatter_genome_size_vs_gc_trna.png")
print(summary_folder / "scatter_genome_size_vs_cds.png")
print(summary_folder / "scatter_genome_size_vs_trna.png")
print(summary_folder / "scatter_genome_size_combined.png")
