#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.io
from pathlib import Path
from PIL import Image
import math
import seaborn as sns

sns.set(style="whitegrid")

# ========= Parámetros =========
TARGET_GSE = "GSE122960"            # Solo este GSE
VIRIDIS = "viridis"                 # Paleta a igualar
MT_VMIN = None                      # Pon aquí el vmin que usaste (ej. 0) o deja None
MT_VMAX = None                      # Pon aquí el vmax que usaste (ej. 20) o deja None

FIG_LIST = [
    "violin_n_genes.png",
    "violin_n_counts.png",
    "violin_pct_counts_mt.png",
    "hist_n_counts.png",
    "hist_n_genes.png",
    "hist_pct_counts_mt.png",
    "scatter_qc.png",
    "top_expressed_genes.png"
]

def mad(x):
    return np.median(np.abs(x - np.median(x)))

# ========= Rutas =========
base = Path(__file__).resolve().parent.parent
raw_base = base / "data" / "processed"
fig_pref_base = base / "figures" / "00.prefiltered"
res_pref_base = base / "results" / "00.prefiltered"  # por si luego quieres guardar algo
fig_pref_base.mkdir(parents=True, exist_ok=True)
res_pref_base.mkdir(parents=True, exist_ok=True)

# ========= Loop sólo GSE122960 =========
gse_path = raw_base / TARGET_GSE
if not gse_path.exists():
    raise FileNotFoundError(f"No encuentro {gse_path}")

print(f"\n📁 PRE-QC para: {TARGET_GSE}")

gse_fig_dir = fig_pref_base / TARGET_GSE
gse_fig_dir.mkdir(parents=True, exist_ok=True)

for gsm in gse_path.iterdir():
    if not gsm.is_dir() or not gsm.name.startswith("GSM"):
        continue

    sample_id = gsm.name
    print(f"\n🔬 Muestra: {sample_id}")

    fig_out = gse_fig_dir / sample_id
    fig_out.mkdir(parents=True, exist_ok=True)

    try:
        # ---- Lectura datos crudos ----
        h5_file = gsm / "matrix.h5"
        if h5_file.exists():
            print("📥 Leyendo H5")
            adata = sc.read_10x_h5(h5_file)
            adata.var_names_make_unique()
        else:
            mtx = gsm / "matrix.mtx.gz"
            genes = gsm / "genes.tsv.gz"
            barcodes = gsm / "barcodes.tsv.gz"
            if not (mtx.exists() and genes.exists() and barcodes.exists()):
                print(f"⚠️  Archivos faltantes en {sample_id}, se omite.")
                continue

            print("📥 Leyendo MTX/TSV")
            matrix = scipy.io.mmread(mtx).tocsc().T
            barcodes_df = pd.read_csv(barcodes, header=None, sep="\t")
            genes_df = pd.read_csv(genes, header=None, sep="\t")

            adata = sc.AnnData(X=matrix)
            adata.obs_names = barcodes_df[0].values
            adata.var_names = genes_df[1] if genes_df.shape[1] > 1 else genes_df[0]
            adata.var_names_make_unique()

        # ---- QC sin filtrar ----
        adata.obs["sample"] = sample_id
        adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1
        adata.obs["n_counts"] = adata.X.sum(axis=1).A1
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        print(f"   🔢 Celdas totales (sin filtrar): {adata.n_obs}")

        # ---- Violinplots ----
        for metric in ["n_genes", "n_counts", "pct_counts_mt"]:
            sc.pl.violin(
                adata,
                keys=[metric],
                jitter=0.4,
                show=False,
                scale='width'
            )
            plt.tight_layout()
            plt.savefig(fig_out / f"violin_{metric}.png", dpi=150)
            plt.close()

        # ---- Histogramas ----
        for metric in ["n_genes", "n_counts", "pct_counts_mt"]:
            plt.figure(figsize=(8, 5))
            plt.hist(adata.obs[metric], bins=50, color="steelblue", edgecolor="black")
            plt.title(f"Histograma de {metric}")
            plt.xlabel(metric)
            plt.ylabel("Nº células")
            plt.tight_layout()
            plt.savefig(fig_out / f"hist_{metric}.png", dpi=150)
            plt.close()

        # ---- Scatter QC con misma paleta ----
        try:
            sc.pl.scatter(
                adata,
                x="n_counts",
                y="n_genes",
                color="pct_counts_mt",
                color_map=VIRIDIS,
                vmin=MT_VMIN,
                vmax=MT_VMAX,
                show=False
            )
            plt.tight_layout()
            plt.savefig(fig_out / "scatter_qc.png", dpi=150)
            plt.close()
        except TypeError:
            # Fallback matplotlib
            plt.figure(figsize=(6, 4))
            cvals = adata.obs["pct_counts_mt"]
            plt.scatter(
                adata.obs["n_counts"],
                adata.obs["n_genes"],
                c=cvals,
                cmap=VIRIDIS,
                vmin=MT_VMIN,
                vmax=MT_VMAX,
                s=3, alpha=0.6, linewidths=0
            )
            plt.xlabel("n_counts"); plt.ylabel("n_genes")
            cb = plt.colorbar(); cb.set_label("pct_counts_mt")
            plt.tight_layout()
            plt.savefig(fig_out / "scatter_qc.png", dpi=150)
            plt.close()

        # ---- Top genes ----
        sc.pl.highest_expr_genes(adata, n_top=20, show=False)
        plt.tight_layout()
        plt.savefig(fig_out / "top_expressed_genes.png", dpi=150)
        plt.close()

    except Exception as e:
        print(f"❌ Error en {sample_id}: {e}")

# ========= Combinar figuras =========
print(f"\n🧩 Combinando figuras para {TARGET_GSE}")

gse_fig_root = gse_fig_dir  # where GSM folders are
combo_out = gse_fig_root / "figures_recopilation"
combo_out.mkdir(exist_ok=True)

# Parte 1: por tipo de figura (grid de GSMs)
for fig_name in FIG_LIST:
    images = []
    titles = []
    for gsm in sorted(os.listdir(gse_fig_root)):
        gsm_path = gse_fig_root / gsm
        if not gsm_path.is_dir() or not gsm.startswith("GSM"):
            continue
        fp = gsm_path / fig_name
        if fp.exists():
            try:
                images.append(Image.open(fp))
                titles.append(gsm)
            except Exception as e:
                print(f"❌ Error abriendo {fp}: {e}")

    if not images:
        print(f"⛔ No hay imágenes para {fig_name}")
        continue

    cols = 4
    rows = math.ceil(len(images) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = axes.flatten()

    for i, (img, title) in enumerate(zip(images, titles)):
        axes[i].imshow(img)
        axes[i].set_title(title, fontsize=9)
        axes[i].axis("off")

    for j in range(len(images), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    out_path = combo_out / f"combined_grid_{fig_name}"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"   ✅ {out_path.name}")

# (Opcional) Barplots resumidos si tienes qc_summary.csv en 00.prefiltered
qc_csv = res_pref_base / TARGET_GSE / "qc_summary.csv"
if qc_csv.exists():
    df = pd.read_csv(qc_csv, index_col=0)
    ncols = 2
    nrows = (len(df.columns) + 1) // 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14, 5 * nrows))
    axes = axes.flatten()

    for i, column in enumerate(df.columns):
        ax = axes[i]
        sns.barplot(x=df.index, y=df[column], ax=ax, color="gray")
        ax.set_title(f"{TARGET_GSE} - {column}", fontsize=12)
        ax.set_ylabel(column)
        ax.set_xlabel("Muestra")
        ax.tick_params(axis='x', rotation=90)
        for bar in ax.patches:
            h = bar.get_height()
            if h > 0:
                ax.annotate(f'{h:.2f}',
                            (bar.get_x() + bar.get_width() / 2, h * 0.05),
                            ha='center', va='bottom',
                            fontsize=8, color='white', rotation=90)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    bar_out = combo_out / f"{TARGET_GSE}_qc_barplots_combined.png"
    plt.savefig(bar_out, dpi=300)
    plt.close()
    print(f"   📊 {bar_out.name}")

print("\n✅ Proceso terminado.")

