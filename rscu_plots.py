#!/usr/bin/env python3

"""
rscu_plots.py

Reads files from ~/Stenophage/RSCU:
- codon_usage_RSCU.csv
- codon_usage_counts.csv
- trna_totals_by_phage.csv

Makes one side-by-side figure:
A. GC12 vs GC3 with regression line
B. RSCU PCA colored by tRNA presence

K279a = green, D457 = red.
Output: ~/Stenophage/RSCU/rscu_gc_pca_combined.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

folder = Path("~/Stenophage/RSCU").expanduser()
rscu_csv = folder / "codon_usage_RSCU.csv"
codon_csv = folder / "codon_usage_counts.csv"
trna_csv = folder / "trna_totals_by_phage.csv"
output_png = folder / "rscu_gc_pca_combined.png"

host_colors = {"K279a": "green", "D457": "red"}

CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]

def load_rscu(path):
    df = pd.read_csv(path).set_index("phage")
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.fillna(df.mean(axis=0))
    return df

def load_trna_presence(path):
    df = pd.read_csv(path)
    df["tRNA_total"] = pd.to_numeric(df["tRNA_total"], errors="coerce").fillna(0)
    return pd.Series(df["tRNA_total"].values > 0, index=df["phage"].astype(str))

def calc_gc12_gc3(row):
    total_12 = 0
    gc_12 = 0
    total_3 = 0
    gc_3 = 0

    for cod in CODONS:
        n = row.get(cod, 0)
        if pd.isna(n):
            n = 0
        n = int(n)

        if n == 0:
            continue

        b1, b2, b3 = cod[0], cod[1], cod[2]

        total_12 += 2 * n
        gc_12 += ((b1 in "GC") + (b2 in "GC")) * n

        total_3 += n
        gc_3 += (b3 in "GC") * n

    if total_12 == 0 or total_3 == 0:
        return pd.Series({"GC12": np.nan, "GC3": np.nan})

    return pd.Series({
        "GC12": gc_12 / total_12,
        "GC3": gc_3 / total_3
    })

def add_fit(ax, x, y):
    m, b = np.polyfit(x, y, 1)
    r = np.corrcoef(x, y)[0, 1]
    r2 = r ** 2
    xs = np.linspace(x.min(), x.max(), 200)

    ax.plot(xs, m * xs + b, color="black", linestyle="--", linewidth=2)

    eq = f"y = {m:.3f}x + {b:.3f}\nR² = {r2:.3f}\nr = {r:.3f}"
    ax.text(0.04, 0.96, eq, transform=ax.transAxes, va="top", ha="left", fontsize=9, bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.85))

def make_gc_plot(ax):
    df = pd.read_csv(codon_csv)
    gc_df = df.apply(calc_gc12_gc3, axis=1)
    df = pd.concat([df[["phage"]], gc_df], axis=1).dropna()

    hosts = df["phage"].isin(["K279a", "D457"])
    phages = ~hosts

    ax.scatter(df.loc[phages, "GC12"], df.loc[phages, "GC3"], s=35, alpha=0.55, color="tab:blue", edgecolors="none", label="Phages")

    for host, color in host_colors.items():
        h = df[df["phage"] == host]
        if len(h) > 0:
            ax.scatter(h["GC12"], h["GC3"], s=80, alpha=0.7, color=color, edgecolors="black", linewidths=0.8, label=host)

    add_fit(ax, df["GC12"].to_numpy(), df["GC3"].to_numpy())

    ax.set_title("A  GC12 vs GC3 content")
    ax.set_xlabel("GC12")
    ax.set_ylabel("GC3")
    ax.legend(frameon=False, fontsize=8)

def make_pca_plot(ax):
    X_df = load_rscu(rscu_csv)
    trna_presence = load_trna_presence(trna_csv)
    has_trna = X_df.index.to_series().astype(str).map(trna_presence).fillna(False).astype(bool)

    X_scaled = StandardScaler().fit_transform(X_df.to_numpy())
    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(X_scaled)

    pc1 = pcs[:, 0]
    pc2 = pcs[:, 1]

    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100

    host_mask = X_df.index.isin(["K279a", "D457"])
    yes_mask = has_trna.to_numpy() & ~host_mask
    no_mask = (~has_trna.to_numpy()) & ~host_mask

    ax.scatter(pc1[no_mask], pc2[no_mask], s=28, alpha=0.6, label="tRNA no", color="tab:orange", edgecolors="none")
    ax.scatter(pc1[yes_mask], pc2[yes_mask], s=28, alpha=0.6, label="tRNA yes", color="tab:blue", edgecolors="none")

    for host, color in host_colors.items():
        if host in X_df.index:
            i = list(X_df.index).index(host)
            ax.scatter(pc1[i], pc2[i], s=80, alpha=0.7, label=host, color=color, edgecolors="black", linewidths=0.8)

    ax.set_title("B  PCA of RSCU profiles")
    ax.set_xlabel(f"PC1 ({var1:.1f}%)")
    ax.set_ylabel(f"PC2 ({var2:.1f}%)")
    ax.legend(frameon=False, fontsize=8)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

make_gc_plot(axes[0])
make_pca_plot(axes[1])

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(output_png, dpi=300)
plt.close(fig)

print("Saved:", output_png)
