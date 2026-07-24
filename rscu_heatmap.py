#!/usr/bin/env python3

"""
rscu_heatmap.py

Reads:
~/Stenophage/RSCU/codon_usage_RSCU.csv

Makes:
- codons on right side
- phages on bottom
- blue → red gradient
- codon labels include amino acid names

Output:
~/Stenophage/RSCU/rscu_heatmap.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

folder = Path("~/Stenophage/RSCU").expanduser()

rscu_csv = folder / "codon_usage_RSCU.csv"
output_png = folder / "rscu_heatmap.png"

AA = {
    "TTT":"Phe","TTC":"Phe",
    "TTA":"Leu","TTG":"Leu",
    "CTT":"Leu","CTC":"Leu",
    "CTA":"Leu","CTG":"Leu",

    "ATT":"Ile","ATC":"Ile","ATA":"Ile",
    "ATG":"Met",

    "GTT":"Val","GTC":"Val",
    "GTA":"Val","GTG":"Val",

    "TCT":"Ser","TCC":"Ser",
    "TCA":"Ser","TCG":"Ser",
    "AGT":"Ser","AGC":"Ser",

    "CCT":"Pro","CCC":"Pro",
    "CCA":"Pro","CCG":"Pro",

    "ACT":"Thr","ACC":"Thr",
    "ACA":"Thr","ACG":"Thr",

    "GCT":"Ala","GCC":"Ala",
    "GCA":"Ala","GCG":"Ala",

    "TAT":"Tyr","TAC":"Tyr",

    "CAT":"His","CAC":"His",
    "CAA":"Gln","CAG":"Gln",

    "AAT":"Asn","AAC":"Asn",
    "AAA":"Lys","AAG":"Lys",

    "GAT":"Asp","GAC":"Asp",
    "GAA":"Glu","GAG":"Glu",

    "TGT":"Cys","TGC":"Cys",
    "TGG":"Trp",

    "CGT":"Arg","CGC":"Arg",
    "CGA":"Arg","CGG":"Arg",
    "AGA":"Arg","AGG":"Arg",

    "GGT":"Gly","GGC":"Gly",
    "GGA":"Gly","GGG":"Gly",

    "TAA":"Stop",
    "TAG":"Stop",
    "TGA":"Stop"
}

df = pd.read_csv(rscu_csv)

df = df.set_index("phage")

df = df.apply(
    pd.to_numeric,
    errors="coerce"
)

df = df.replace(
    [np.inf, -np.inf],
    np.nan
)

df = df.dropna(
    axis=1,
    how="all"
)

df = df.fillna(
    df.mean(axis=0)
)

plot_df = df.T

# order codons by average RSCU
mean_rscu = plot_df.mean(axis=1)

plot_df = plot_df.loc[
    mean_rscu.sort_values(
        ascending=False
    ).index
]

new_labels = []

for cod in plot_df.index:

    aa = AA.get(
        cod,
        "?"
    )

    new_labels.append(
        f"{cod} ({aa})"
    )

plot_df.index = new_labels

plt.figure(
    figsize=(
        max(
            12,
            len(plot_df.columns)*0.22
        ),
        16
    )
)

ax = sns.heatmap(
    plot_df,
    cmap="coolwarm",
    vmin=0,
    vmax=2,
    xticklabels=True,
    yticklabels=True,
    cbar_kws={
        "label":"RSCU"
    }
)

ax.yaxis.tick_right()
ax.yaxis.set_label_position(
    "right"
)

ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title(
    "RSCU Heatmap"
)

plt.xticks(
    rotation=90,
    fontsize=7
)

plt.yticks(
    rotation=0,
    fontsize=8
)

plt.tight_layout()

plt.savefig(
    output_png,
    dpi=300
)

plt.close()

print(
    "Saved:",
    output_png
)
