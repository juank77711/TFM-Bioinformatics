import scanpy as sc
import numpy as np
from sklearn.metrics import adjusted_rand_score
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# ========= Configuración =========
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
input_file = os.path.join(base_dir, "results", "04.integration", "all_samples_umap_50k.h5ad")
results_dir = os.path.join(base_dir, "results", "05.clustering")
figures_dir = os.path.join(base_dir, "figures", "05.clustering")
os.makedirs(results_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

resolutions = [0.3, 0.5]
n_iter = 10
sample_size = 5000
random_seed = 42

# ========= Cargar datos =========
print(f"📥 Cargando archivo: {input_file}")
adata_full = sc.read_h5ad(input_file)
np.random.seed(random_seed)

ari_results = {r: [] for r in resolutions}

# ========= ARI: Estabilidad de resolución =========
for i in range(n_iter):
    print(f"🔁 Iteración {i+1}/{n_iter}")
    sampled_cells = np.random.choice(adata_full.obs_names, size=sample_size, replace=False)
    adata_sub = adata_full[sampled_cells].copy()
    sc.pp.neighbors(adata_sub, use_rep="X_pca_harmony")

    cluster_labels = {}

    for res in resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata_sub, resolution=res, key_added=key)
        cluster_labels[res] = adata_sub.obs[key].values.copy()

        if i > 0:
            prev_labels = prev_labels_dict[res]
            common_idx = np.intersect1d(sampled_cells, prev_cells)
            if len(common_idx) > 0:
                ari = adjusted_rand_score(
                    cluster_labels[res][np.isin(sampled_cells, common_idx)],
                    prev_labels[np.isin(prev_cells, common_idx)]
                )
                ari_results[res].append(ari)

    prev_labels_dict = cluster_labels
    prev_cells = sampled_cells

# ========= Guardar resultados ARI =========
ari_df = pd.DataFrame(dict([(f"res_{r}", pd.Series(vals)) for r, vals in ari_results.items()]))
ari_csv_path = os.path.join(results_dir, "ari_scores.csv")
ari_df.to_csv(ari_csv_path)
print(f"📝 Resultados ARI guardados en: {ari_csv_path}")

# ========= Plot ARI =========
ari_df_melted = ari_df.melt(var_name="Resolution", value_name="ARI")
plt.figure(figsize=(8, 5))
sns.boxplot(x="Resolution", y="ARI", data=ari_df_melted)
plt.title("Stability (ARI) across resolutions")
plt.ylabel("Adjusted Rand Index (ARI)")
plt.xlabel("Leiden resolution")
plt.tight_layout()
plot_path = os.path.join(figures_dir, "ari_resolution_comparison.png")
plt.savefig(plot_path, dpi=150)
plt.close()
print(f"📊 Gráfico ARI guardado en: {plot_path}")

# ========= Selección de resolución óptima =========
median_aris = {r: np.median(vals) for r, vals in ari_results.items()}
best_res = max(median_aris, key=median_aris.get)
print(f"\n✅ Resolución seleccionada automáticamente: {best_res} (mediana ARI = {median_aris[best_res]:.4f})")

# ========= Clustering final completo =========
leiden_key = "leiden"
print("📡 Calculando vecinos para clustering final...")
sc.pp.neighbors(adata_full, use_rep='X_pca_harmony')
sc.tl.leiden(adata_full, resolution=best_res, key_added=leiden_key)

# ========= UMAP final por clúster =========
print("🖼️ Generando UMAP por clúster final...")
sc.pl.umap(adata_full, color=leiden_key, title=f"Leiden Clusters (res={best_res})", show=False)
umap_path = os.path.join(figures_dir, f"umap_leiden_final_res_{str(best_res).replace('.', '_')}.png")
plt.savefig(umap_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"🖼️ UMAP final guardado en: {umap_path}")

# ========= Guardar objeto AnnData =========
output_file = os.path.join(results_dir, "all_samples_clustered_50k.h5ad")
adata_full.write(output_file)
print(f"💾 AnnData con clustering final guardado en: {output_file}")
print("🎉 Proceso completo con éxito.")

