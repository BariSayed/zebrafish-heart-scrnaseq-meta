# Run from the repository root; all paths below are repository-relative.
# ============================================================================
# Figure 3A | Integrated UMAP overview, 6 panels (cluster/cell-type/GSE/dpi/condition + legend)
# INPUT:  umap_full_export.csv.gz + new_cluster_annotation.csv (cluster->cell_type)
# OUTPUT: Figure_3A.png  (assembled with Figure3B into manuscript Figure 3)
# NOTE: df + colour maps (originally defined in upstream cells) are inlined here so the script is standalone.
# ============================================================================
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
try:
    from figure_style import apply_figure_style
except Exception:
    def apply_figure_style():  # no-op fallback if the figure-style helper is unavailable
        pass
# ---- load UMAP export + annotate cell_type ----
df  = pd.read_csv("03_processed_data/clustering/umap_full_export.csv.gz", compression='gzip')
ann = pd.read_csv("03_processed_data/clustering/new_cluster_annotation.csv")
df['cell_type'] = df['cluster'].map(dict(zip(ann['new_cluster'], ann['final_cell_type'])))
assert df['cell_type'].isna().sum()==0, "unmapped clusters"
# ---- cluster colour map (tab20 + tab20b, 26) ----
clusters_sorted = sorted(df['cluster'].unique()); n = len(clusters_sorted)
palette = np.vstack([plt.cm.tab20(np.linspace(0,1,20)),
                     plt.cm.tab20b(np.linspace(0,1,20))[[2,6,10,14,18,19]]])[:n]
cluster_color = {cl: palette[i] for i, cl in enumerate(clusters_sorted)}
# ---- cell-type colour map (ordered by abundance; tab20 + Set3) ----
celltype_order = sorted(df['cell_type'].unique(), key=lambda x: -df[df['cell_type']==x].shape[0])
ct_palette = np.vstack([plt.cm.tab20(np.linspace(0,1,20)), plt.cm.Set3(np.linspace(0,1,12))])[:len(celltype_order)]
celltype_color = {ct: ct_palette[i] for i, ct in enumerate(celltype_order)}

apply_figure_style()
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, wspace=0.05, hspace=0.25)

def scatter_panel(ax, color_col, color_map, title, legend=True, legend_ncol=1, s=0.8, sort_groups=None):
    groups = sort_groups if sort_groups else sorted(df[color_col].unique())
    for g in groups:
        sub = df[df[color_col]==g]
        ax.scatter(sub['UMAP_1'], sub['UMAP_2'], s=s, color=color_map[g], linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, fontsize=11, loc='left')
    return groups

# a) cluster ID
ax = fig.add_subplot(gs[0,0])
scatter_panel(ax, 'cluster', cluster_color, "a   Cluster ID", s=0.6)
for cl in clusters_sorted:
    sub = df[df['cluster']==cl]
    ax.text(sub['UMAP_1'].median(), sub['UMAP_2'].median(), str(cl), fontsize=6.5, ha='center', va='center',
            fontweight='bold', path_effects=[])

# b) cell type (collapsed coloring, no legend here - legend shown separately for space)
ax = fig.add_subplot(gs[0,1])
scatter_panel(ax, 'cell_type', celltype_color, "b   Cell type", s=0.6)

# c) GSE / dataset
gse_colors = {"GSE227190":"#8fd3c6","GSE159032":"#fdf5b0","GSE262689":"#c5b8db"}
ax = fig.add_subplot(gs[0,2])
scatter_panel(ax, 'GSE', gse_colors, "c   Source dataset (GSE)", s=0.6)
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=gse_colors[g], markersize=7, label=g) for g in gse_colors]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.0,1.0), frameon=False, fontsize=8)

# d) dpi
dpi_colors = {"0dpi":"#1f77b4","1dpi":"#ff7f0e","3dpi":"#2ca02c","7dpi":"#d62728","14dpi":"#9467bd"}
ax = fig.add_subplot(gs[1,0])
scatter_panel(ax, 'dpi', dpi_colors, "d   Days post-injury", s=0.6, sort_groups=["0dpi","1dpi","3dpi","7dpi","14dpi"])
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=dpi_colors[g], markersize=7, label=g) for g in ["0dpi","1dpi","3dpi","7dpi","14dpi"]]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.0,1.0), frameon=False, fontsize=8)

# e) condition
cond_colors = {"Healthy":"#2ca02c","Injured":"#d62728"}
ax = fig.add_subplot(gs[1,1])
scatter_panel(ax, 'Condition', cond_colors, "e   Condition", s=0.6)
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=cond_colors[g], markersize=7, label=g) for g in cond_colors]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.0,1.0), frameon=False, fontsize=8)

# f) cell-type legend panel (text-only, since 23 colors won't fit as scatter legend cleanly)
ax = fig.add_subplot(gs[1,2]); ax.axis('off')
ax.set_title("f   Cell type legend", fontsize=11, loc='left')
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=celltype_color[ct], markersize=7, label=ct) for ct in celltype_order]
ax.legend(handles=handles, loc='center left', frameon=False, fontsize=6.5, labelspacing=0.6, ncol=1)

fig.savefig("05_figures/main/Figure_3A.png", dpi=300, bbox_inches='tight')
plt.close(fig)
print("saved")
