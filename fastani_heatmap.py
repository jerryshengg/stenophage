import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

out_dir = os.path.expanduser("~/Stenophage/ANI")
input_file = os.path.join(out_dir, "fastani.csv")

df = pd.read_csv(
    input_file,
    sep=r"\s+",
    engine="python",
    header=None
)

df.columns = ["query", "reference", "ani", "fragments_aligned", "total_fragments"]

df["query"] = df["query"].str.replace("_query", "", regex=False)
df["reference"] = df["reference"].str.replace("_ref", "", regex=False)

ani_matrix = df.pivot(index="query", columns="reference", values="ani")

all_phages = sorted(set(ani_matrix.index).union(set(ani_matrix.columns)))
ani_matrix = ani_matrix.reindex(index=all_phages, columns=all_phages)

np.fill_diagonal(ani_matrix.values, 100)

# ensure white figure background and seaborn style
sns.set(style="white")
plt.rcParams["figure.facecolor"] = "white"

# make symmetric for clustering
ani_symmetric = ani_matrix.combine_first(ani_matrix.T)
ani_symmetric = ani_symmetric.combine_first(ani_symmetric.T)

# missing values dark/low
ani_filled = ani_symmetric.fillna(70)

g = sns.clustermap(
    ani_filled,
    cmap="YlGnBu",  # yellow -> green -> blue
    vmin=75,
    vmax=100,
    figsize=(22, 22),
    linewidths=0.1,
    linecolor="gray",
    method="average",
    metric="euclidean",
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "ANI (%)"}
)

# make sure the figure and heatmap axes use a white background (not dark colormap end)
g.fig.patch.set_facecolor("white")
g.ax_heatmap.set_facecolor("white")

# labels like old heatmap
g.ax_heatmap.set_title("Clustered Heatmap of Stenotrophomonas Phages (Full Figure)", pad=80)
g.ax_heatmap.set_xlabel("Reference genome")
g.ax_heatmap.set_ylabel("Query genome")

g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    rotation=90,
    fontsize=5
)
g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_yticklabels(),
    rotation=0,
    fontsize=5
)

# Output directory inside the user's home Stenophage folder
out_dir = os.path.expanduser("~/Stenophage/ANI")
os.makedirs(out_dir, exist_ok=True)

# Save high-resolution PNG (lossless) at high DPI; JPEG outputs removed per request
out_base = os.path.join(out_dir, "fastani_clustered_heatmap_labeled")
png_path = out_base + ".png"
g.fig.savefig(png_path, dpi=600, bbox_inches="tight")
print(f"Saved high-resolution PNG: {png_path}")

plt.show()

# ----------------------
# Second, zoomed heatmap (top-left block)
# ----------------------
# Use the clustered, reordered matrix produced by the clustermap
data2d = g.data2d

# Assumption: the 'interesting' high-ANI region is grouped in the top-left after clustering.
# We'll take a square top-left block. Default block size chosen to produce fewer squares for
# a simplified/zoomed figure. Change `block_size` as needed.
block_size = 43
block_size = min(block_size, data2d.shape[0], data2d.shape[1])

zoom_df = data2d.iloc[:block_size, :block_size]

# Create a new figure for the zoomed heatmap
fig2, ax2 = plt.subplots(figsize=(10, 10))
sns.heatmap(
    zoom_df,
    ax=ax2,
    cmap="YlGnBu",
    vmin=75,
    vmax=100,
    linewidths=0.1,
    linecolor="gray",
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "ANI (%)"},
)

ax2.set_title(f"Clustered Heatmap of High-ANI Stenotrophomonas Phages", pad=16)
ax2.set_xlabel("Reference genome")
ax2.set_ylabel("Query genome")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, fontsize=6)
ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0, fontsize=6)

# Save the zoomed outputs at high resolution (PNG only)
out_base_zoom = os.path.join(out_dir, "fastani_clustered_heatmap_labeled_zoom")
png_zoom = out_base_zoom + ".png"

fig2.savefig(png_zoom, dpi=600, bbox_inches="tight")
print(f"Saved zoomed PNG: {png_zoom}")

plt.show()
