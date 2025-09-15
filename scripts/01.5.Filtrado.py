import os
from PIL import Image
import matplotlib.pyplot as plt
import math
import pandas as pd
import seaborn as sns

sns.set(style="whitegrid")

# Figuras a combinar por GSM
fig_names = [
    "violin_n_genes.png",
    "violin_n_counts.png",
    "violin_pct_counts_mt.png",
    "hist_n_counts.png",
    "hist_n_genes.png",
    "hist_pct_counts_mt.png",
    "scatter_qc.png",
    "top_expressed_genes.png"
]

# Ruta base del proyecto
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
figures_root = os.path.join(project_root, "figures", "01.filtered")
results_root = os.path.join(project_root, "results", "01.filtered")

# Recorrer cada GSE
for gse in os.listdir(figures_root):
    gse_path = os.path.join(figures_root, gse)
    if not os.path.isdir(gse_path) or not gse.startswith("GSE"):
        continue

    print(f"\n📁 Procesando GSE: {gse}")
    output_dir = os.path.join(gse_path, "figures_recopilation")
    os.makedirs(output_dir, exist_ok=True)

    # Parte 1: Combinar imágenes por tipo
    for fig_name in fig_names:
        images = []
        titles = []

        for gsm in os.listdir(gse_path):
            gsm_path = os.path.join(gse_path, gsm)
            if not os.path.isdir(gsm_path) or not gsm.startswith("GSM"):
                continue

            fig_path = os.path.join(gsm_path, fig_name)
            if os.path.exists(fig_path):
                try:
                    img = Image.open(fig_path)
                    images.append(img)
                    titles.append(gsm)
                except Exception as e:
                    print(f"❌ Error al abrir {fig_path}: {e}")

        if not images:
            print(f"⛔ No se encontraron imágenes para {fig_name}, se omite.")
            continue

        cols = 4
        rows = math.ceil(len(images) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
        axes = axes.flatten()

        for i, (img, title) in enumerate(zip(images, titles)):
            axes[i].imshow(img)
            axes[i].set_title(title, fontsize=10)
            axes[i].axis("off")

        for j in range(len(images), len(axes)):
            axes[j].axis("off")

        plt.tight_layout()
        output_path = os.path.join(output_dir, f"combined_grid_{fig_name}")
        plt.savefig(output_path, dpi=150)
        plt.close()
        print(f"✅ Figura combinada guardada: {output_path}")

    # Parte 2: Barplots de resumen desde qc_summary.csv
    qc_path = os.path.join(results_root, gse, "qc_summary.csv")
    if not os.path.exists(qc_path):
        print(f"⚠️  No se encontró qc_summary.csv en results/{gse}, se omite resumen QC.")
        continue

    df = pd.read_csv(qc_path, index_col=0)
    num_cols = len(df.columns)
    ncols = 2
    nrows = (num_cols + 1) // 2

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14, 5 * nrows))
    axes = axes.flatten()

    for i, column in enumerate(df.columns):
        ax = axes[i]
        sns.barplot(x=df.index, y=df[column], ax=ax, color='gray')
        ax.set_title(f"{gse} - {column}", fontsize=12)
        ax.set_ylabel(column)
        ax.set_xlabel("Muestra")
        ax.tick_params(axis='x', rotation=90)

        for bar in ax.patches:
            height = bar.get_height()
            if height > 0:
                ax.annotate(
                    f'{height:.2f}',
                    (bar.get_x() + bar.get_width() / 2, height * 0.05),
                    ha='center', va='bottom',
                    fontsize=8, color='white', rotation=90
                )

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    qc_output_path = os.path.join(output_dir, f"{gse}_qc_barplots_combined.png")
    plt.savefig(qc_output_path, dpi=300)
    plt.close()
    print(f"📊 Barplot resumen QC guardado: {qc_output_path}")

