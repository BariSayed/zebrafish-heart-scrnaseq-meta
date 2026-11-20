# Run from the repository root; all paths below are repository-relative.
# ============================================================================
# Figure 8 | Pseudotime distribution boxplots by days-post-injury, 5 panels
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  5 *_pseudotime_values.csv (from Figure7A/7B)
# OUTPUT:    Figure_8.png
# ============================================================================
import pandas as pd, numpy as np, matplotlib.pyplot as plt
panels=[
 ("A","Cardiomyocytes","03_processed_data/pseudotime/Cardiomyocytes_pseudotime_values.csv"),
 ("B","Endocardial cell","03_processed_data/pseudotime/Endocardial_cell_pseudotime_values.csv"),
 ("C","Epicardial cell","03_processed_data/pseudotime/Epicardial_cell_pseudotime_values.csv"),
 ("D","Fibroblasts","03_processed_data/pseudotime/Fibroblasts_pseudotime_values.csv"),
 ("E","Macrophages","03_processed_data/pseudotime/Macrophages_pseudotime_values.csv"),
]
order=["0dpi","1dpi","3dpi","7dpi","14dpi"]

fig,axes=plt.subplots(2,3,figsize=(15,9)); axes=axes.ravel()
for i,(L,name,vid) in enumerate(panels):
    ax=axes[i]
    pt=pd.read_csv(vid); pt=pt[np.isfinite(pt['pseudotime'])]
    data=[pt[pt.dpi==d]['pseudotime'].values for d in order]; ns=[len(x) for x in data]
    bp=ax.boxplot(data,tick_labels=[d.replace('dpi','') for d in order],showfliers=False,patch_artist=True,widths=0.6)
    for patch in bp['boxes']: patch.set_facecolor('#8fb4d9'); patch.set_edgecolor('#2b4a6f')
    for med in bp['medians']: med.set_color('#c0392b'); med.set_linewidth(1.6)
    ax.set_title(name,fontsize=12); ax.set_xlabel('days post-injury',fontsize=10); ax.set_ylabel('pseudotime',fontsize=10)
    ax.text(-0.08,1.04,L,transform=ax.transAxes,fontsize=18,fontweight='bold',va='top')
    ax.spines[['top','right']].set_visible(False)
    ymax=ax.get_ylim()[1]
    for j,n in enumerate(ns): ax.text(j+1,ymax*0.98,f"n={n}",ha='center',va='top',fontsize=7,color='#555')
axes[5].axis('off'); fig.tight_layout()
fig.savefig("05_figures/main/Figure_8.png",dpi=300,bbox_inches='tight')
print("saved")
