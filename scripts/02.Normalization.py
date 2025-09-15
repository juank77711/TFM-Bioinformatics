import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from PIL import Image, ImageOps
import math

# Establecer rutas
base_dir = Path(__file__).resolve().parent.parent
filtered_base = base_dir / "results" / "01.filtered"
results_base = base_dir / "results" / "02.normalized"
figures_base = base_dir / "figures" / "02.normalized"

# Iterar sobre los GSE
for gse_dir in filtered_base.iterdir():
    if not gse_dir.is_dir() or not gse_dir.name.startswith("GSE"):
        continue

    print(f"\n🔍 Procesando {gse_dir.name}...")

    # Crear carpetas
    gse_results = results_base / gse_dir.name
    gse_figures = figures_base / gse_dir.name
    gse_results.mkdir(parents=True, exist_ok=True)
    gse_figures.mkdir(parents=True, exist_ok=True)

    # CSV resumen
    summary_file = gse_results / "qc_normalized_summary.csv"
    if not summary_file.exists():
        pd.DataFrame(columns=["SampleID", "n_cells", "mean_genes", "mean_counts"]).to_csv(summary_file, index=False)

    # Track imágenes exitosas y fallidas
    images_successful = []
    titles_successful = []
    images_failed = []
    titles_failed = []

    # Iterar sobre GSM
    for gsm_dir in gse_dir.iterdir():
        if not gsm_dir.is_dir() or not gsm_dir.name.startswith("GSM"):
            continue

        sample_id = gsm_dir.name
        data_file = gsm_dir / "filtered_data.h5ad"
        if not data_file.exists():
            print(f"⚠️  No se encontró filtered_data.h5ad en {sample_id}")
            images_failed.append(None)
            titles_failed.append(sample_id)
            continue

        try:
            print(f"🔄 Normalizando muestra: {sample_id}")
            adata = sc.read_h5ad(data_file)

            if adata.shape[0] == 0 or adata.shape[1] == 0:
                print(f"⚠️  Muestra vacía: {sample_id}, saltando.")
                images_failed.append(None)
                titles_failed.append(sample_id)
                continue

            # Normalización
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat')

            # Rutas de salida
            result_out = gse_results / sample_id
            figure_out = gse_figures / sample_id
            result_out.mkdir(parents=True, exist_ok=True)
            figure_out.mkdir(parents=True, exist_ok=True)

            # Guardar figura HVG
            fig_path = figure_out / "highly_variable_genes_plot.png"
            sc.pl.highly_variable_genes(adata, show=False)
            plt.savefig(fig_path)
            plt.close()

            # Guardar AnnData
            adata.write(result_out / "normalized_data.h5ad")

            # CSV resumen
            row = {
                "SampleID": sample_id,
                "n_cells": adata.n_obs,
                "mean_genes": adata.obs["n_genes"].mean(),
                "mean_counts": adata.obs["n_counts"].mean(),
            }
            pd.DataFrame([row]).to_csv(summary_file, mode='a', header=False, index=False)

            # Guardar para combinar
            images_successful.append(Image.open(fig_path))
            titles_successful.append(sample_id)

        except Exception as e:
            print(f"❌ Error en {sample_id}: {e}")
            images_failed.append(None)
            titles_failed.append(sample_id)
            continue

    # ========== COMBINAR FIGURAS ==========
    if images_successful or images_failed:
        all_titles = titles_successful + titles_failed
        all_images = images_successful + [ImageOps.grayscale(Image.new("RGB", images_successful[0].size)) if img is None else img for img in images_failed]

        n_imgs = len(all_images)
        cols = 4
        rows = math.ceil(n_imgs / cols)

        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
        axes = axes.flatten()

        for i, (img, title) in enumerate(zip(all_images, all_titles)):
            axes[i].imshow(img)
            axes[i].set_title(title, fontsize=10)
            axes[i].axis("off")

        # Desactivar ejes vacíos
        for j in range(len(all_images), len(axes)):
            axes[j].axis("off")

        plt.tight_layout()

        # Carpeta para figuras combinadas
        recop_path = gse_figures / "figures_recopilation"
        recop_path.mkdir(parents=True, exist_ok=True)
        save_path = recop_path / "combined_hvg.png"
        plt.savefig(save_path, dpi=150)
        plt.close()

        print(f"📊 Figura recopilada guardada en: {save_path}")

print("\n✅ Normalización + figura combinada completada.")

