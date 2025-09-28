import os
import scanpy as sc
import harmonypy as hm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_hex

# ====================
# 📁 Rutas y ajustes
# ====================
base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
input_file = os.path.join(base, "results", "03.concatenated", "new_all_samples_concatenated_normalized_ready.h5ad")
results_dir = os.path.join(base, "results", "04.integration")
figures_dir = os.path.join(base, "figures", "04.integration")
os.makedirs(results_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

corrected_file = os.path.join(results_dir, "all_samples_corrected_qc.h5ad")
output_file    = os.path.join(results_dir, "all_samples_integrated_comparison.h5ad")
umap_output    = os.path.join(results_dir, "all_samples_umap_50k.h5ad")
umap_fig_base  = os.path.join(figures_dir, "umap")

# Puntos más “intensos”
POINT_SIZE_MAIN = 14     # comparativa pre/post
POINT_ALPHA_MAIN = 0.95
POINT_SIZE_INDV = 12     # UMAPs individuales
POINT_ALPHA_INDV = 0.95

sc.settings.set_figure_params(dpi=120, frameon=False)

# === helpers para marco negro y guardado ===
def add_black_frame(ax):
    ax.set_frame_on(True)  # Fuerza visibilidad del marco
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.4)

def save_fig_with_black_frames(fig, axes, outpath, dpi=220):
    if isinstance(axes, (list, np.ndarray)):
        for ax in np.ravel(axes):
            add_black_frame(ax)
    else:
        add_black_frame(axes)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

# ====================
# 📥 Lectura inicial
# ====================
print(f"📥 Cargando archivo: {input_file}")
adata = sc.read_h5ad(input_file)
adata.obs_names_make_unique()

# Asegurar campos categóricos
adata.obs['sample'] = adata.obs.get('sample', adata.obs.get('sample_id'))
adata.obs['sample'] = adata.obs['sample'].astype("category")
adata.obs['GSM']    = adata.obs['GSM'].astype("category")
adata.obs['GSE']    = adata.obs['GSE'].astype("category")

# ========================
# 💾 Guardar QC corregido
# ========================
adata.write(corrected_file)
print(f"💾 Archivo corregido guardado: {corrected_file}")

# ===============
# 🔄 Re-cargar
# ===============
adata = sc.read_h5ad(corrected_file)

# ===============
# ✨ Submuestreo
# ===============
max_cells = 50000
if adata.n_obs > max_cells:
    print(f"⚠️ Submuestreando a {max_cells} células...")
    np.random.seed(42)
    adata = adata[np.random.choice(adata.obs_names, size=max_cells, replace=False)].copy()

# ============================
# 🎨 Paletas fijas (GSM y GSE)
# ============================
gsm_list = list(adata.obs['GSM'].cat.categories)
gse_list = list(adata.obs['GSE'].cat.categories)

palette_gsm = sns.color_palette("husl", len(gsm_list))
palette_gse = sns.color_palette("tab20", len(gse_list))

adata.uns['GSM_colors'] = [to_hex(c) for c in palette_gsm]
adata.uns['GSE_colors'] = [to_hex(c) for c in palette_gse]

# ===============
# 📊 PCA & Harmony (por GSM)
# ===============
print("🔍 HVG + PCA (batch_key=GSM)…")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", batch_key="GSM")
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
adata.obsm['X_pca_original'] = adata.obsm['X_pca'].copy()

print("🌀 UMAP pre-Harmony…")
sc.pp.neighbors(adata, use_rep='X_pca_original', key_added='neighbors_pre')
sc.tl.umap(adata, neighbors_key='neighbors_pre')
adata.obsm['X_umap_pre'] = adata.obsm['X_umap'].copy()

print("🔗 Harmony por GSM…")
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs[['GSM']], vars_use=['GSM'])
adata.obsm['X_pca_harmony'] = ho.Z_corr.T

print("🌀 UMAP post-Harmony…")
sc.pp.neighbors(adata, use_rep='X_pca_harmony', key_added='neighbors_post')
sc.tl.umap(adata, neighbors_key='neighbors_post')
adata.obsm['X_umap_post'] = adata.obsm['X_umap'].copy()

# ================================
# 🧠 Diccionario de condiciones
# ================================
condition_map = {
    'GSM4827199':'Sano','GSM4850115':'COVID-19','GSM4827188':'IPF','GSM3891626':'IPF',
    'GSM4698183':'COVID-19','GSM4698178':'COVID-19','GSM4850119':'COVID-19','GSM4827180':'Ssc-ILD',
    'GSM4850116':'COVID-19','GSM3891620':'Sano','GSM3489198':'Ssc-ILD','GSM4850112':'COVID-19',
    'GSM3891623':'Sano','GSM4698177':'COVID-19','GSM4827205':'IPF','GSM4827173':'Ssc-ILD',
    'GSM3891628':'IPF','GSM4850120':'COVID-19','GSM3489195':'Sano','GSM3489194':'Ssc-ILD',
    'GSM3489185':'Sano','GSM3891629':'IPF','GSM4827189':'IPF','GSM3660648':'Sano',
    'GSM3489197':'Sano','GSM3666100':'Sano','GSM3489191':'Sano','GSM4516282':'COVID-19',
    'GSM4827172':'Ssc-ILD','GSM4516281':'COVID-19','GSM4827185':'Sano','GSM4827169':'Ssc-ILD',
    'GSM4827193':'IPF','GSM4850111':'COVID-19','GSM3891633':'Ssc-ILD','GSM4827179':'Ssc-ILD',
    'GSM3489193':'Sano','GSM4850117':'COVID-19','GSM4827164':'Sano','GSM3489196':'Ssc-ILD',
    'GSM4827182':'Sano','GSM4827191':'IPF','GSM3660655':'IPF','GSM4827178':'Ssc-ILD',
    'GSM4827195':'IPF','GSM4827168':'Ssc-ILD','GSM4827176':'Ssc-ILD','GSM4827194':'IPF',
    'GSM3660654':'IPF','GSM3489188':'IPF','GSM4827174':'Ssc-ILD','GSM3660656':'IPF',
    'GSM4516280':'COVID-19','GSM4850114':'COVID-19','GSM3666099':'Sano','GSM4698182':'COVID-19',
    'GSM4827198':'Sano','GSM4827201':'IPF','GSM3891621':'Sano','GSM4698184':'COVID-19',
    'GSM4827187':'IPF','GSM3666097':'Sano','GSM3660645':'Sano','GSM4827196':'Sano',
    'GSM3660657':'IPF','GSM4850118':'COVID-19','GSM4827161':'IPF','GSM3489187':'Sano',
    'GSM4827159':'IPF','GSM4698179':'COVID-19','GSM4827184':'Sano','GSM4516279':'COVID-19',
    'GSM3666104':'Ssc-ILD','GSM3660646':'Sano','GSM4827190':'IPF','GSM4827200':'Sano',
    'GSM3666103':'Ssc-ILD','GSM3660651':'IPF','GSM4827202':'IPF','GSM3666102':'Ssc-ILD',
    'GSM3891632':'Ssc-ILD','GSM4850113':'COVID-19','GSM3891622':'Sano','GSM3666105':'Ssc-ILD',
    'GSM3660642':'Sano','GSM3660653':'IPF','GSM3666098':'Sano','GSM4827166':'Sano',
    'GSM4827160':'IPF','GSM4827170':'Ssc-ILD','GSM4827197':'Sano','GSM4827171':'Ssc-ILD',
    'GSM3489182':'Sano','GSM4827162':'IPF','GSM3660652':'IPF','GSM4827183':'Sano',
    'GSM4827203':'IPF','GSM4827204':'IPF','GSM3660658':'IPF','GSM3666106':'Ssc-ILD',
    'GSM4827192':'IPF','GSM4827167':'Sano','GSM3660647':'Sano','GSM4827175':'Ssc-ILD',
    'GSM4698181':'COVID-19','GSM3489189':'Sano','GSM3666096':'Sano','GSM3660644':'Sano',
    'GSM3660641':'Sano','GSM3891627':'IPF','GSM3489190':'IPF','GSM3489184':'IPF',
    'GSM3666101':'Ssc-ILD','GSM4698180':'COVID-19','GSM3666108':'Ssc-ILD','GSM4827177':'Ssc-ILD',
    'GSM4698176':'COVID-19','GSM4827165':'Sano','GSM4827163':'Sano','GSM3489183':'IPF',
    'GSM4827158':'IPF','GSM4827181':'Sano','GSM3660649':'Sano','GSM3666107':'Ssc-ILD',
    'GSM3660643':'Sano','GSM3660650':'Sano','GSM4827186':'IPF'
}
adata.obs['condition'] = adata.obs['sample_id'].map(condition_map).fillna("Desconocido")
adata.obs['condition'] = adata.obs['condition'].astype("category")

# =====================
# 💾 Guardar objeto
# =====================
adata.write(output_file)
print(f"💾 Objeto final integrado guardado: {output_file}")

# ========================
# 📸 Comparativa Pre/Post coloreada por GSE
# ========================
print("📸 UMAP Pre vs Post Harmony (coloreado por GSE)…")
fig, axes = plt.subplots(1, 2, figsize=(22, 9))
sc.pl.embedding(
    adata, basis='X_umap_pre', color='GSE', ax=axes[0], show=False,
    title='Pre-Harmony by GSE', legend_loc='right margin', frameon=False,
    size=POINT_SIZE_MAIN, alpha=POINT_ALPHA_MAIN
)
sc.pl.embedding(
    adata, basis='X_umap_post', color='GSE', ax=axes[1], show=False,
    title='Post-Harmony by GSE', legend_loc='right margin', frameon=False,
    size=POINT_SIZE_MAIN, alpha=POINT_ALPHA_MAIN
)
save_fig_with_black_frames(fig, axes, f"{umap_fig_base}_pre_post_harmony_by_GSE_colored_fixed.png")

# ========================
# ✨ Submuestreo final
# ========================
print("🔬 Submuestreo final para visualización (50k)…")
adata_sub = adata.copy()
adata_sub.write(umap_output)
print(f"💾 Submuestreo guardado: {umap_output}")

# ========================
# 🎨 UMAPs individuales
# ========================
print("📊 UMAPs individuales…")
for key, fname in [
    ("condition",   "umap_by_condition.png"),
    ("n_genes",     "umap_by_n_genes.png"),
    ("pct_counts_mt","umap_by_mito.png"),
]:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata_sub, color=key, ax=ax, show=False,
        legend_loc='right margin', frameon=False,
        size=POINT_SIZE_INDV, alpha=POINT_ALPHA_INDV
    )
    save_fig_with_black_frames(fig, ax, os.path.join(figures_dir, fname))

print("\n✅ Integración por GSM con bordes negros en todos los PNGs — listo.")
