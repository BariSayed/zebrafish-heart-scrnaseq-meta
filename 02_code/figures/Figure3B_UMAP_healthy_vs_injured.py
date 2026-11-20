# Run from the repository root; all paths below are repository-relative.
# ============================================================================
# Figure 3B | UMAP split Healthy vs Injured, 8 major cell types (others greyed)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  umap_full_export.csv
# OUTPUT:    Figure_3B.png
# NOTE: Lower section of the final Figure 4 (assembled with Figure3A).
# ============================================================================
# Use the 8 major interpretable cell types as in old Figure 4 (not all 26 clusters - too busy for a split panel)
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
# load UMAP export + annotate cell_type so this script runs standalone (same as Figure3A)
df  = pd.read_csv("03_processed_data/clustering/umap_full_export.csv.gz", compression="gzip")
ann = pd.read_csv("03_processed_data/clustering/new_cluster_annotation.csv")
df["cell_type"] = df["cluster"].map(dict(zip(ann["new_cluster"], ann["final_cell_type"])))

major_map = {
    "Endocardial cell":"Endocardial cell", "Endocardial cell (frzb+)":"Endocardial cell",
    "Cardiomyocytes":"Cardiomyocytes", "Cardiomyocytes (ventricular/sarcomeric)":"Cardiomyocytes",
    "Epicardial cells":"Epicardial cells", "Epicardial cell":"Epicardial cells",
    "Macrophages":"Macrophages",
    "Fibroblasts":"Fibroblasts",
    "T cells / lymphocytes":"T cells",
    "Neutrophils":"Neutrophils",
    "Perivascular cell / pericyte":"Perivascular cell",
}
df['major_type'] = df['cell_type'].map(major_map)  # NaN for others (erythroid/progenitor/etc - greyed out, matches old fig scope)

major_order = ["Cardiomyocytes","Endocardial cell","Epicardial cells","Macrophages","Fibroblasts","T cells","Neutrophils","Perivascular cell"]
major_colors = {ct: plt.cm.tab10(i) for i,ct in enumerate(major_order)}

fig, axes = plt.subplots(1,2, figsize=(14,6.5), sharex=True, sharey=True)
for ax, cond in zip(axes, ["Healthy","Injured"]):
    sub_all = df[df['Condition']==cond]
    # background: non-major cells in light grey
    bg = sub_all[sub_all['major_type'].isna()]
    ax.scatter(bg['UMAP_1'], bg['UMAP_2'], s=0.4, color='#e8e8e8', linewidths=0, rasterized=True, zorder=1)
    for ct in major_order:
        s2 = sub_all[sub_all['major_type']==ct]
        ax.scatter(s2['UMAP_1'], s2['UMAP_2'], s=0.8, color=major_colors[ct], linewidths=0, rasterized=True, zorder=2, label=ct)
    ax.set_title(cond, fontsize=13)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_xlabel("UMAP 1", fontsize=10); 
axes[0].set_ylabel("UMAP 2", fontsize=10)
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=major_colors[ct], markersize=9, label=ct) for ct in major_order]
handles.append(plt.Line2D([0],[0], marker='o', color='w', markerfacecolor='#cccccc', markersize=9, label='Other (erythroid/progenitor/etc.)'))
axes[1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02,1.0), frameon=False, fontsize=9)
fig.suptitle("Major cell-type distribution: Healthy vs Injured", fontsize=13, y=1.02)
fig.subplots_adjust(wspace=0.05)
fig.savefig("05_figures/main/Figure_3B.png", dpi=300, bbox_inches='tight')
plt.close(fig)
print("saved; n major-typed cells:", df['major_type'].notna().sum(), "/ 77375")
