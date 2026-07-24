#!/usr/bin/env python3

"""
plot_phage_summary.py

Reads ~/Stenophage/Summary/phage_summary.csv.
Creates GC cutoff plots and p-value stats.
Outputs everything to ~/Stenophage/Summary.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy import stats

summary_folder = Path("~/Stenophage/Summary").expanduser()
csv_file = summary_folder / "phage_summary.csv"
outdir = summary_folder
gc_cutoff = 59.0

NUMERIC_COLS = ["genome_length_bp", "GC_percent", "CDS_count", "tRNA_count", "CDS_per_10kb"]
LOW_COLOR = "tab:orange"
HIGH_COLOR = "tab:blue"

def thousands(x, pos):
    try:
        return f"{int(x):,}"
    except Exception:
        return str(x)

def load_df(csv_path):
    df = pd.read_csv(csv_path)
    for c in NUMERIC_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def require_columns(df, cols, csv_path):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in {csv_path}: {missing}")

def split_by_gc_cutoff(df, gc_cutoff, cols):
    needed = list(set(["GC_percent", *cols]))
    d = df[needed].dropna().copy()
    low = d[d["GC_percent"] < gc_cutoff].copy()
    high = d[d["GC_percent"] >= gc_cutoff].copy()
    return low, high

def two_group_tests(df, value_col, gc_cutoff, do_mannwhitney=True):
    low, high = split_by_gc_cutoff(df, gc_cutoff, cols=[value_col])
    x = low[value_col].to_numpy(dtype=float)
    y = high[value_col].to_numpy(dtype=float)

    out = {
        "metric": value_col,
        "gc_cutoff": gc_cutoff,
        "n_low": int(x.size),
        "n_high": int(y.size),
        "mean_low": float(np.nanmean(x)) if x.size else np.nan,
        "mean_high": float(np.nanmean(y)) if y.size else np.nan,
        "welch_t_p": np.nan,
        "mannwhitney_p": np.nan
    }

    if x.size >= 2 and y.size >= 2:
        _, p = stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
        out["welch_t_p"] = float(p)

        if do_mannwhitney:
            try:
                _, p2 = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            except TypeError:
                _, p2 = stats.mannwhitneyu(x, y, alternative="two-sided")
            out["mannwhitney_p"] = float(p2)

    return out

def scatter_gc_vs_trna_cutoff(ax, df, gc_cutoff, title):
    low, high = split_by_gc_cutoff(df, gc_cutoff, cols=["tRNA_count"])

    ax.scatter(low["GC_percent"], low["tRNA_count"], s=26, alpha=0.9, color=LOW_COLOR, label=f"GC < {gc_cutoff:g}")
    ax.scatter(high["GC_percent"], high["tRNA_count"], s=26, alpha=0.9, color=HIGH_COLOR, label=f"GC ≥ {gc_cutoff:g}")

    ax.axvline(x=gc_cutoff, linestyle=":", linewidth=2)
    ax.set_title(title)
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("tRNA count")
    ax.legend(frameon=False, fontsize=9, loc="best")

def boxplot_by_cutoff_with_points(ax, df, value_col, gc_cutoff, title, ylabel, format_thousands=False, point_size=14, jitter=0.08):
    low, high = split_by_gc_cutoff(df, gc_cutoff, cols=[value_col])
    low_vals = low[value_col].to_numpy()
    high_vals = high[value_col].to_numpy()

    x_low, x_high = 1.0, 2.0
    rng = np.random.default_rng(0)

    if low_vals.size:
        ax.scatter(rng.normal(x_low, jitter, size=low_vals.size), low_vals, s=point_size, alpha=0.55, color=LOW_COLOR, zorder=1)
    if high_vals.size:
        ax.scatter(rng.normal(x_high, jitter, size=high_vals.size), high_vals, s=point_size, alpha=0.55, color=HIGH_COLOR, zorder=1)

    bp = ax.boxplot([low_vals, high_vals], positions=[x_low, x_high], widths=0.55, showfliers=False, patch_artist=True, zorder=2)

    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor([LOW_COLOR, HIGH_COLOR][i])
        patch.set_alpha(0.35)
        patch.set_linewidth(1.2)

    for line in bp["whiskers"] + bp["caps"]:
        line.set_linewidth(1.2)

    for i, med in enumerate(bp["medians"]):
        med.set_linewidth(1.6)
        med.set_color([LOW_COLOR, HIGH_COLOR][i])

    ax.set_xticks([x_low, x_high])
    ax.set_xticklabels([f"<{gc_cutoff:g}", f"≥{gc_cutoff:g}"])
    ax.set_title(title)
    ax.set_ylabel(ylabel)

    if format_thousands:
        ax.yaxis.set_major_formatter(FuncFormatter(thousands))

def make_scatter_gc_vs_trna(df, outpath, gc_cutoff):
    fig, ax = plt.subplots(1, 1, figsize=(6.8, 5.2))
    scatter_gc_vs_trna_cutoff(ax, df, gc_cutoff, title="GC content (%) vs tRNA count")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

def make_boxplots_by_gc_cutoff_1x3(df, outpath, gc_cutoff):
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), constrained_layout=True)

    boxplot_by_cutoff_with_points(axes[0], df, "genome_length_bp", gc_cutoff, "A  Genome length", "bp", format_thousands=True)
    boxplot_by_cutoff_with_points(axes[1], df, "CDS_count", gc_cutoff, "B  ORF count", "ORFs")
    boxplot_by_cutoff_with_points(axes[2], df, "CDS_per_10kb", gc_cutoff, "C  CDS per 10 kb", "CDS / 10 kb")

    fig.savefig(outpath, dpi=300)
    plt.close(fig)

def main():
    df = load_df(csv_file)

    require_columns(df, ["GC_percent", "tRNA_count"], csv_file)
    require_columns(df, ["genome_length_bp", "CDS_count", "CDS_per_10kb"], csv_file)

    fig_scatter = outdir / "figure_scatter_gc_vs_trna.png"
    fig_boxes = outdir / "figure_boxplots_by_gc_cutoff.png"
    stats_csv = outdir / "gc_stats.csv"

    make_scatter_gc_vs_trna(df, fig_scatter, gc_cutoff)
    make_boxplots_by_gc_cutoff_1x3(df, fig_boxes, gc_cutoff)

    metrics = ["genome_length_bp", "CDS_count", "CDS_per_10kb"]
    rows = [two_group_tests(df, m, gc_cutoff, do_mannwhitney=True) for m in metrics]
    pd.DataFrame(rows).to_csv(stats_csv, index=False)

    print(f"Saved scatter plot: {fig_scatter}")
    print(f"Saved boxplots: {fig_boxes}")
    print(f"Saved stats: {stats_csv}")

if __name__ == "__main__":
    main()
