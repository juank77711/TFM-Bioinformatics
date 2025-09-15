import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.io
from pathlib import Path

def mad(x):
    return np.median(np.abs(x - np.median(x)))

# Rutas base
base = Path(__file__).resolve().parent.parent
raw_base = base / "data" / "processed"
figures_base = base / "figures" / "01.filtered"
results_base = base / "results" / "01.filtered"

for gse in raw_base.iterdir():
    if not gse.is_dir() or not gse.name.startswith("GSE"):
        continue

    print(f"\n📁 Procesando GSE: {gse.name}")
    gse_figures = figures_base / gse.name
    gse_results = results_base / gse.name
    gse_figures.mkdir(parents=True, exist_ok=True)
    gse_results.mkdir(parents=True, exist_ok=True)

    summary_file = gse_results / "qc_summary.csv"
    if not summary_file.exists():
        pd.DataFrame(columns=["SampleID", "n_cells_final", "mean_genes", "mean_counts", "mean_pct_mt"]).to_csv(summary_file, index=False)

    for gsm in gse.iterdir():
        if not gsm.is_dir() or not gsm.name.startswith("GSM"):
            continue

        print(f"\n🔬 Procesando muestra: {gsm.name}")
        sample_id = gsm.name
        fig_out = gse_figures / sample_id
        res_out = gse_results / sample_id
        fig_out.mkdir(parents=True, exist_ok=True)
        res_out.mkdir(parents=True, exist_ok=True)

        try:
            # Buscar archivos estandarizados
            h5_file = gsm / "matrix.h5"
            adata = None

            if h5_file.exists():
                print(f"📥 Leyendo archivo H5: matrix.h5")
                adata = sc.read_10x_h5(h5_file)
                adata.var_names_make_unique()
            else:
                mtx = gsm / "matrix.mtx.gz"
                genes = gsm / "genes.tsv.gz"
                barcodes = gsm / "barcodes.tsv.gz"

                print(f"📥 Leyendo archivos MTX/TSV: {mtx.name}, {genes.name}, {barcodes.name}")
                if not (mtx.exists() and genes.exists() and barcodes.exists()):
                    print(f"⚠️  Archivos faltantes en {gsm.name}, se omite.")
                    continue

                matrix = scipy.io.mmread(mtx).tocsc().T
                barcodes_df = pd.read_csv(barcodes, header=None, sep="\t")
                genes_df = pd.read_csv(genes, header=None, sep="\t")

                adata = sc.AnnData(X=matrix)
                adata.obs_names = barcodes_df[0].values
                adata.var_names = genes_df[1] if genes_df.shape[1] > 1 else genes_df[0]
                adata.var_names_make_unique()

            # QC
            adata.obs["sample"] = sample_id
            adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1
            adata.obs["n_counts"] = adata.X.sum(axis=1).A1
            adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

            print(f"   🔢 Celdas antes de filtrar: {adata.n_obs}")

            # Umbrales
            med_genes = np.median(adata.obs["n_genes"])
            mad_genes = mad(adata.obs["n_genes"])
            upper_genes = med_genes + 3 * mad_genes if mad_genes > 0 else 2500
            lower_genes = 200

            if np.isnan(adata.obs["pct_counts_mt"]).all():
                mt_thresh = 20
            else:
                med_mt = np.median(adata.obs["pct_counts_mt"])
                mad_mt = mad(adata.obs["pct_counts_mt"])
                mt_thresh = med_mt + 3 * mad_mt if mad_mt > 0 else 20

            print(f"   ➤ Umbrales: n_genes [{lower_genes}-{upper_genes}], pct_mt < {mt_thresh:.2f}")

            # Filtrado
            adata = adata[adata.obs["n_genes"] > lower_genes, :]
            adata = adata[adata.obs["n_genes"] < upper_genes, :]
            adata = adata[adata.obs["n_counts"] > 0, :]
            adata = adata[adata.obs["pct_counts_mt"] < mt_thresh, :]
            sc.pp.filter_genes(adata, min_cells=3)

            print(f"   ✅ Celdas tras filtrar: {adata.n_obs}")
            if adata.n_obs == 0:
                print("   ⚠️  Ninguna célula pasó el filtrado.")
                continue

            # Guardar AnnData y resumen
            adata.write(res_out / "filtered_data.h5ad")

            summary_row = {
                "SampleID": sample_id,
                "n_cells_final": adata.n_obs,
                "mean_genes": adata.obs["n_genes"].mean(),
                "mean_counts": adata.obs["n_counts"].mean(),
                "mean_pct_mt": adata.obs["pct_counts_mt"].mean()
            }
            pd.DataFrame([summary_row]).to_csv(summary_file, mode='a', header=False, index=False)

            # Gráficos
            for metric in ["n_genes", "n_counts", "pct_counts_mt"]:
                sc.pl.violin(adata, [metric], jitter=0.4, show=False, scale='width')
                plt.tight_layout()
                plt.savefig(fig_out / f"violin_{metric}.png")
                plt.close()

                plt.figure(figsize=(8, 5))
                plt.hist(adata.obs[metric], bins=50, color="steelblue", edgecolor="black")
                plt.title(f"Histograma de {metric}")
                plt.xlabel(metric)
                plt.ylabel("Nº células")
                plt.tight_layout()
                plt.savefig(fig_out / f"hist_{metric}.png")
                plt.close()

            sc.pl.scatter(adata, x="n_counts", y="n_genes", color="pct_counts_mt", show=False)
            plt.tight_layout()
            plt.savefig(fig_out / "scatter_qc.png")
            plt.close()

            sc.pl.highest_expr_genes(adata, n_top=20, show=False)
            plt.tight_layout()
            plt.savefig(fig_out / "top_expressed_genes.png")
            plt.close()

            print(f"   📦 Guardado: {res_out}/filtered_data.h5ad")

        except Exception as e:
            print(f"❌ Error en {gsm.name}: {e}")

