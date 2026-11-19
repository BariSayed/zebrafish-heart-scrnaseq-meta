# ============================================================================
# Figure 6B | Curated GO/KEGG dot-plot, 4 panels (GO-up, KEGG-up, GO-down, KEGG-down)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  GO_BP/KEGG_enrichment_{UP,DOWN}_per_cluster.csv (from Figure6A)
# OUTPUT:    Figure_6.png
# ============================================================================
# Figure 8 (GO/KEGG dot-plot, 4 panels). PART 1 = data load + curated terms; PART 2 = panel assembly.

# ===== PART 1: data + curated term lists =====
import pandas as pd, numpy as np, matplotlib.pyplot as plt
apply_figure_style()

# reload corrected CSVs (cluster-12 fixed versions)
go_up = pd.read_csv(host.artifact_path("964c8697-203d-4ec4-9688-763b91fc707c"))
go_dn = pd.read_csv(host.artifact_path("8840618c-7c21-41f9-bc36-50da60243816"))
kg_up = pd.read_csv(host.artifact_path("84e6e8fe-8f78-4085-8e02-e45986118cdd"))
kg_dn = pd.read_csv(host.artifact_path("28661ff2-d1d2-46f3-ba77-ad06a076e297"))

col_order = [("2","Cardiomyocyte"),("5","Cardiomyocyte (vent.)"),("0","Endocardial"),("8","Endocardial"),
    ("3","Epicardial"),("22","Epicardial"),("19","Fibroblast"),("21","Perivascular"),("14","Vascular endo."),
    ("6","Macrophage"),("18","Neutrophil"),("10","T cell"),("13","Erythroid")]

def best_match(df, cl_id, substrs):
    sub = df[df['Cluster'].str.startswith(cl_id+'_')]
    mask = sub['Description'].apply(lambda d: any(s.lower() in d.lower() for s in substrs))
    hit = sub[mask]
    if len(hit)==0: return None
    r = hit.nsmallest(1,'p.adjust').iloc[0]
    return (-np.log10(r['p.adjust']), int(r['Count']))

# curated, split strictly by source database
GO_UP = [("Extracellular matrix organization",["extracellular matrix organization"]),
         ("Skeletal system development",["skeletal system development"]),
         ("Actin filament polymerization",["actin filament polymerization","actin filament capping","barbed-end actin"]),
         ("Endothelial cell migration",["endothelial cell migration"]),
         ("Erythrocyte differentiation",["erythrocyte differentiation","erythrocyte homeostasis"])]
KEGG_UP = [("ECM-receptor interaction",["ECM-receptor interaction"]),
           ("Focal adhesion / integrin",["Focal adhesion","Integrin signaling"]),
           ("Protein processing in ER",["Protein processing in endoplasmic reticulum"]),
           ("Glycolysis / carbon metabolism",["Glycolysis / Gluconeogenesis","Carbon metabolism","Biosynthesis of amino acids"]),
           ("Xenobiotic / glutathione metabolism",["Glutathione metabolism","cytochrome P450","Metabolism of xenobiotics"])]
GO_DN = [("RNA splicing",["RNA splicing"]),
         ("Protein folding (chaperone)",["protein folding","protein refolding","protein maturation"]),
         ("Response to heat / metal ion",["response to heat","response to temperature","response to metal ion","response to cadmium"])]
KEGG_DN = [("Spliceosome",["Spliceosome"]),
           ("Protein processing in ER",["Protein processing in endoplasmic reticulum"]),
           ("Integrated stress response",["Integrated stress response"])]

def build(df, groups):
    rows=[]
    for label,subs in groups:
        for cl_id,ct in col_order:
            m=best_match(df,cl_id,subs)
            if m: rows.append(dict(pathway=label,col=f"{ct} ({cl_id})",neglog=m[0],count=m[1]))
    return pd.DataFrame(rows)

d_gu,d_ku,d_gd,d_kd = build(go_up,GO_UP),build(kg_up,KEGG_UP),build(go_dn,GO_DN),build(kg_dn,KEGG_DN)
print("GO-up",len(d_gu),"KEGG-up",len(d_ku),"GO-down",len(d_gd),"KEGG-down",len(d_kd))


# ===== PART 2: 4-panel assembly =====
col_seq=[f"{ct} ({c})" for c,ct in col_order]
groups_lab={"GO_UP":[g[0] for g in GO_UP],"KEGG_UP":[g[0] for g in KEGG_UP],
            "GO_DN":[g[0] for g in GO_DN],"KEGG_DN":[g[0] for g in KEGG_DN]}

def panel(ax,df,ordered,title,cmap):
    cols=[c for c in col_seq if c in set(df['col'])]
    xp={c:i for i,c in enumerate(cols)}; yp={p:i for i,p in enumerate(reversed(ordered))}
    xs=[];ys=[];sz=[];cl=[]
    for _,r in df.iterrows():
        xs.append(xp[r['col']]);ys.append(yp[r['pathway']]);sz.append(r['count']);cl.append(r['neglog'])
    sc=ax.scatter(xs,ys,s=[s*20 for s in sz],c=cl,cmap=cmap,edgecolors='black',linewidths=0.4,zorder=3,vmin=0)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,rotation=45,ha='right',fontsize=7.5)
    ax.set_yticks(range(len(ordered))); ax.set_yticklabels(list(reversed(ordered)),fontsize=8)
    ax.set_xlim(-0.6,len(cols)-0.4); ax.set_ylim(-0.6,len(ordered)-0.4)
    ax.grid(color='#ececec',linewidth=0.5,zorder=0)
    ax.set_title(title,fontsize=11,loc='left',pad=6)
    for sp in ['top','right']: ax.spines[sp].set_visible(False)
    cb=plt.colorbar(sc,ax=ax,fraction=0.04,pad=0.02); cb.set_label("-log10 adj. p",fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

fig,axes=plt.subplots(2,2,figsize=(14,8.5))
panel(axes[0,0],d_gu,groups_lab["GO_UP"],"a   GO BP — upregulated","YlOrRd")
panel(axes[0,1],d_ku,groups_lab["KEGG_UP"],"b   KEGG — upregulated","YlOrRd")
panel(axes[1,0],d_gd,groups_lab["GO_DN"],"c   GO BP — downregulated","PuBu")
panel(axes[1,1],d_kd,groups_lab["KEGG_DN"],"d   KEGG — downregulated","PuBu")
for cv in [3,10,20]:
    axes[0,1].scatter([],[],s=cv*20,c='grey',edgecolors='black',linewidths=0.4,label=str(cv))
axes[0,1].legend(title="genes",loc='upper left',bbox_to_anchor=(1.12,1.0),frameon=False,fontsize=7,title_fontsize=7.5,labelspacing=1.1)
fig.subplots_adjust(left=0.13,right=0.93,top=0.95,bottom=0.17,wspace=0.75,hspace=0.6)
fig.savefig("Figure_6.png",dpi=300,bbox_inches='tight')
plt.close(fig)
print("saved")
