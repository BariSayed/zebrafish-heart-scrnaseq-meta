# Run from the repository root; all paths below are repository-relative.
# ============================================================================
# Figure 7C | Pseudotime on UMAP, 5-panel composite (plots the 7A+7B outputs)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  umap coords CSV + 5 *_pseudotime_values.csv (from 6A/6B)
# OUTPUT:    Figure_7.png
# ============================================================================
import pandas as pd, numpy as np, matplotlib.pyplot as plt
um=pd.read_csv("03_processed_data/clustering/umap_full_export.csv.gz", usecols=['cell','UMAP_1','UMAP_2'])
panels=[
 ("A","Cardiomyocytes","03_processed_data/pseudotime/Cardiomyocytes_pseudotime_values.csv"),
 ("B","Endocardial cell","03_processed_data/pseudotime/Endocardial_cell_pseudotime_values.csv"),
 ("C","Epicardial cell","03_processed_data/pseudotime/Epicardial_cell_pseudotime_values.csv"),
 ("D","Fibroblasts","03_processed_data/pseudotime/Fibroblasts_pseudotime_values.csv"),
 ("E","Macrophages","03_processed_data/pseudotime/Macrophages_pseudotime_values.csv"),
]
fig,axes=plt.subplots(2,3,figsize=(15,9))
axes=axes.ravel()
for i,(L,name,vid) in enumerate(panels):
    ax=axes[i]
    pt=pd.read_csv(vid)
    m=pt.merge(um,on='cell',how='inner')
    m=m[np.isfinite(m['pseudotime'])]
    sc=ax.scatter(m['UMAP_1'],m['UMAP_2'],c=m['pseudotime'],cmap='plasma',s=3,linewidths=0)
    fig.colorbar(sc,ax=ax,shrink=0.75,label='pseudotime')
    ax.set_title(f"{name} (n={len(m):,})",fontsize=12,loc='center')
    ax.set_xlabel('UMAP 1',fontsize=10); ax.set_ylabel('UMAP 2',fontsize=10)
    ax.text(-0.08,1.04,L,transform=ax.transAxes,fontsize=18,fontweight='bold',va='top')
    ax.spines[['top','right']].set_visible(False)
axes[5].axis('off')
fig.tight_layout()
fig.savefig("05_figures/main/Figure_7.png",dpi=300,bbox_inches='tight')
print("saved Figure_7 (pseudotime UMAP, 5 panels)")
