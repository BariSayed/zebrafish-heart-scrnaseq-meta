# ============================================================================
# Figure 7C | Pseudotime on UMAP, 5-panel composite (plots the 7A+7B outputs)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  umap coords CSV + 5 *_pseudotime_values.csv (from 6A/6B)
# OUTPUT:    Figure_7.png
# ============================================================================
import pandas as pd, numpy as np, matplotlib.pyplot as plt
um=pd.read_csv(host.artifact_path("db5e8775-175b-4179-879e-aace40fa9250"), usecols=['cell','UMAP_1','UMAP_2'])
panels=[
 ("A","Cardiomyocytes","8a460d04-2cf0-4cdd-9b97-e813cb3bd49b"),
 ("B","Endocardial cell","24f5a792-d229-4d2e-816c-4e3c679ba60a"),
 ("C","Epicardial cell","f0d6edc3-edd6-4b4e-ae83-21187b06fa95"),
 ("D","Fibroblasts","00a562df-32f3-4f69-a201-279f23111efd"),
 ("E","Macrophages","898978fe-b85c-4cdc-b44d-ec619363f02a"),
]
fig,axes=plt.subplots(2,3,figsize=(15,9))
axes=axes.ravel()
for i,(L,name,vid) in enumerate(panels):
    ax=axes[i]
    pt=pd.read_csv(host.artifact_path(vid))
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
fig.savefig("Figure_7.png",dpi=300,bbox_inches='tight')
print("saved Figure_7 (pseudotime UMAP, 5 panels)")
