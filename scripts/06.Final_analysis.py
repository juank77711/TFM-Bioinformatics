#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import celltypist

warnings.filterwarnings("ignore", category=UserWarning)
sc.settings.verbosity = 2

# ============================ Paths ============================
base_dir   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
input_file = os.path.join(base_dir, "results", "05.clustering", "all_samples_clustered_50k.h5ad")
out_dir    = os.path.join(base_dir, "results", "06.final_analysis")
fig_dir    = os.path.join(base_dir, "figures", "06.final_analysis")
os.makedirs(out_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)
sc.settings.figdir = fig_dir

print(f"📥 Cargando: {input_file}")
adata = sc.read_h5ad(input_file)

# ===================== Anotación por condición =================
cond_col = "condition"
if cond_col not in adata.obs.columns:
    adata.obs[cond_col] = "Desconocido"

condition_to_model = {
    "Sano":        "Human_Lung_Atlas.pkl",
    "COVID-19":    "Autopsy_COVID19_Lung.pkl",
    "IPF":         "Human_IPF_Lung.pkl",
    "Ssc-ILD":     "Human_PF_Lung.pkl",
    "Desconocido": "Human_Lung_Atlas.pkl",
}

all_labels = pd.Series(index=adata.obs_names, dtype="object")
for cond, model_name in condition_to_model.items():
    idx = (adata.obs[cond_col].astype(str) == cond).values
    if idx.sum() == 0:
        continue
    print(f"🔍 Anotando {cond} con {model_name} ({idx.sum()} células)")
    model = celltypist.models.Model.load(model=model_name)
    preds = celltypist.annotate(adata[idx].copy(), model=model, majority_voting=True)
    lab = preds.predicted_labels["majority_voting"]
    all_labels.loc[adata.obs_names[idx]] = lab.values

adata.obs["cell_type_final"] = all_labels.astype("category")

# =================== Unificar etiquetas duplicadas ===================
unify_map = {
    "B": "B cells", "B Cells": "B cells", "b cells": "B cells",
    "Mast": "Mast cells", "mast": "Mast cells",
    "Macrophage": "Macrophages", "macrophage": "Macrophages",
    "Macrophages": "Macrophages", "macrophages": "Macrophages",
    "myofibroblast": "Myofibroblasts", "myofibroblasts": "Myofibroblasts",
    "Pericyte": "Pericytes", "pericytes": "Pericytes",
    "SMC": "Smooth muscle cells", "smooth muscle": "Smooth muscle cells",
    "Plasma cells": "Plasma cells", "Plasmablasts": "Plasma cells", "B_Plasma": "Plasma cells",
    "AT1": "AT1", "AT2": "AT2", "ATII": "AT2", "atii": "AT2", "at2": "AT2",
    "basal": "Basal", "Basal resting": "Basal",
}
adata.obs["cell_type_unified"] = (
    adata.obs["cell_type_final"].astype(str).map(lambda x: unify_map.get(x, x))
).astype("category")

# =================== Filtrado de tipos con pocas células ===================
min_cells = 50
counts = adata.obs["cell_type_unified"].value_counts()
keep_types = counts[counts >= min_cells].index.tolist()
print(f"✅ Manteniendo {len(keep_types)} tipos ≥{min_cells} células, resto → 'Other'")

cell_type_filtered = adata.obs["cell_type_unified"].astype(str)
cell_type_filtered = np.where(cell_type_filtered.isin(keep_types), cell_type_filtered, "Other")
adata.obs["cell_type_filtered"] = pd.Categorical(cell_type_filtered)

# =================== UMAP por tipo filtrado + por condición ===================
if "X_umap" not in adata.obsm_keys():
    if "neighbors" not in adata.uns_keys():
        sc.pp.neighbors(adata, use_rep="X_pca_harmony" if "X_pca_harmony" in adata.obsm_keys() else None)
    sc.tl.umap(adata)

# 👇 Leyenda a la derecha (no encima del UMAP)
sc.pl.umap(
    adata, color="cell_type_filtered", legend_loc="right margin",
    legend_fontsize=6, frameon=False, show=False, save="_umap_filtered_unified.png"
)

sc.pl.umap(
    adata, color=cond_col, legend_loc="right margin", frameon=False,
    show=False, save="_umap_by_condition.png"
)

# =================== RGG (genes marcadores por Leiden) ===================
if "leiden" not in adata.obs.columns:
    sc.tl.leiden(adata, resolution=0.5)

need_rgg = True
if "rank_genes_groups" in adata.uns:
    if adata.uns["rank_genes_groups"].get("params", {}).get("groupby", None) == "leiden":
        need_rgg = False

if need_rgg:
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon",
        n_genes=50, use_raw=False
    )

# ============== HEATMAP TOP3 por Leiden (custom, genes en Y) ==============
print("🎨 Heatmap top3 por Leiden con separación vertical…")
names = adata.uns["rank_genes_groups"]["names"]
groups = names.dtype.names
top_genes_per_cluster = {g: list(names[g][:3]) for g in groups}
genes_list = []
for g in groups:
    genes_list.extend([x for x in top_genes_per_cluster[g] if x in adata.var_names])
seen = set()
genes_unique = [x for x in genes_list if not (x in seen or seen.add(x))]

expr = adata.to_df()[genes_unique]
expr["leiden"] = adata.obs["leiden"].values
mean_by_cluster = expr.groupby("leiden").mean().T  # genes x clusters

mean_by_cluster = (mean_by_cluster - mean_by_cluster.mean(axis=1)[:, None]) / (
    mean_by_cluster.std(axis=1)[:, None] + 1e-9
)

col_order = sorted(mean_by_cluster.columns, key=lambda x: int(x))
mean_by_cluster = mean_by_cluster[col_order]

plt.figure(figsize=(max(10, len(col_order) * 0.6), max(8, len(genes_unique) * 0.35)))
ax = sns.heatmap(
    mean_by_cluster,
    cmap="viridis", center=0,
    cbar_kws={"label": "Z-score (media por cluster)"},
    linewidths=0.3, linecolor="black"
)
ax.set_xlabel("leiden", fontsize=12)
ax.set_ylabel("")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)  # separación visual
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "heatmap_top3_by_leiden_unified.png"), dpi=300)
plt.close()

# =================== DOTPLOT TOP10 de RGG (filtrado) ===================
print("🎯 Dotplot top10 (RGG)…")
top10 = []
for g in groups:
    top10.extend([x for x in list(names[g][:10]) if x in adata.var_names])
top10 = list(dict.fromkeys(top10))

sc.pl.dotplot(
    adata, var_names=top10, groupby="cell_type_filtered",
    standard_scale="var", swap_axes=True,
    show=False, save="_dotplot_top10_filtered_unified.png"
)

# =================== DOTPLOT DE CONDICIÓN (genes Y, condición X) ===================
print("🧪 Dotplot por tipo-condición (genes Y, condición X)…")
adata.obs["type_condition"] = (
    adata.obs["cell_type_filtered"].astype(str) + " | " + adata.obs[cond_col].astype(str)
)

marker_genes = top10[:70] if len(top10) >= 70 else top10

groups_tc = pd.Series(adata.obs["type_condition"].unique())
order = groups_tc.sort_values(key=lambda s: s.str.split(" \| ").str[1]).tolist()

sc.pl.dotplot(
    adata,
    var_names=marker_genes,
    groupby="type_condition",
    categories_order=order,
    standard_scale="var",
    swap_axes=True,   # genes en Y, condición en X (vía type|condition)
    show=False,
    save="_dotplot_markers_by_type_condition_swapped.png"
)

# === Heatmap de proporciones en UNA sola figura (sin abreviar nombres) ===
print("📊 Heatmap de proporciones (figura única, etiquetas legibles)…")

# Matriz de proporciones por condición y tipo
prop = (
    adata.obs.groupby(["condition", "cell_type_filtered"]).size()
      .unstack(fill_value=0)
)
prop = prop.div(prop.sum(axis=1), axis=0).round(2)

import textwrap
def wrap_labels(labels, width=12):
    return ['\n'.join(textwrap.wrap(str(l), width=width, break_long_words=False))
            for l in labels]

n_cols = prop.shape[1]
fig_w = max(28, 0.5 * n_cols)   # aumenta 0.5 -> 0.6 si aún hay solapes
fig_h = 7

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

vmax = max(0.05, prop.values.max())

hm = sns.heatmap(
    prop, cmap="Blues", vmin=0, vmax=vmax,
    annot=True, fmt=".1f", annot_kws={"size":8},
    linewidths=0.6, linecolor="lightgray",
    cbar=True, cbar_kws={"label": "Proporción", "shrink": 0.9},
    square=False, ax=ax
)

ax.set_xlabel("cell_type_filtered")
ax.set_ylabel("condition")

# Y: horizontal, tamaño legible y sin solapes
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right", fontsize=12)
ax.tick_params(axis='y', pad=8)

# --- Etiquetas del eje X: verticales, sin dividir, sin wrap ---
n_cols = prop.shape[1]
fig_w = max(32, 0.7 * n_cols)   # más ancho para que quepan
fig.set_size_inches(fig_w, fig_h)

# Reemplaza guiones bajos por espacios para que se lean mejor
xlabels = [str(c).replace('_', ' ') for c in prop.columns]

ax.set_xticklabels(
    xlabels,
    rotation=90,          # vertical
    ha="center",          # centradas bajo cada columna
    va="top",
    fontsize=8            # baja a 7 si aún se tocan
)

# Más aire abajo para las etiquetas largas
plt.subplots_adjust(left=0.16, right=0.98, bottom=0.42, top=0.95)

out_png = os.path.join(fig_dir, "heatmap_modified_unified_SINGLE.png")
plt.savefig(out_png, dpi=300)
plt.close()

print(f"✅ Guardado: {out_png}")

# =================== Guardar objeto ===================
out_h5ad = os.path.join(out_dir, "all_samples_annotated_celltype_unified.h5ad")
adata.write(out_h5ad)
print(f"\n✅ Guardado: {out_h5ad}")

