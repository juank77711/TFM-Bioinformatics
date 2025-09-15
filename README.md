# 🫁 Análisis integrativo de scRNA-seq pulmonar (TFM)

Repositorio asociado al Trabajo Fin de Máster: **“Integración multiestudio de datos de scRNA-seq para la generación de un atlas celular pulmonar: un enfoque técnico-metodológico reproducible”**.

> **Autor:** Juan Carlos Aguado Seminario · **Directora:** Sheila Zúñiga Trejos  
> **Titulación:** Máster en Bioinformática (2024–2025) · **Fecha:** 1 Septiembre 2025

---

## Contenido

### 1. `00.Ordenar.sh`
- **Función**:  
  Organización inicial de carpetas y estandarización de nombres/ubicación de archivos descargados (GEO).
- **Incluye**:
  - Creación de jerarquía `data/`, `results/` y `figures/`.
  - Movido de matrices (`matrix.h5` o `matrix.mtx.gz + genes.tsv.gz + barcodes.tsv.gz`) a `data/processed/GSE*/GSM*/`.
- **Notas**:  
  Edita rutas internas si tu estructura local difiere.

---

### 2. `Prefiltered.py`
- **Función**:  
  **Pre-QC sin filtrado** para un GSE específico (vista rápida de calidad antes del pipeline).
- **Incluye**:
  - Lectura automática de 10x (`*.h5`) o MTX/TSV.
  - Cálculo de métricas: `n_genes`, `n_counts`, `pct_counts_mt`.
  - Figuras: violines, histogramas, scatter QC y genes más expresados por GSM.
  - **Grids combinados** de figuras por tipo.
- **Notas**:  
  Cambia `TARGET_GSE = "GSE122960"` si quieres otra serie. Salidas en `figures/00.prefiltered/`.

---

### 3. `01.Filtrado.py`
- **Función**:  
  **QC + filtrado adaptativo por muestra (GSM)**.
- **Incluye**:
  - Lectura estandarizada (10x H5 o MTX/TSV).
  - Cálculo de métricas QC y marcaje mitocondrial (`MT-*`).
  - **Umbrales**:
    - `n_genes < 200` → descarta.
    - `n_genes > mediana + 3×MAD` → descarta (posibles dobletes).
    - `pct_counts_mt > mediana + 3×MAD` (fallback 20% si MAD=0) → descarta.
    - `n_counts == 0` → descarta.
    - Filtrado de genes con `min_cells = 3`.
  - **Salidas**:
    - `results/01.filtered/<GSE>/<GSM>/filtered_data.h5ad`
    - `figures/01.filtered/<GSE>/<GSM>/*` (violines, histogramas, scatter, top genes)
    - `results/01.filtered/<GSE>/qc_summary.csv` (acumulativo)
- **Notas**:  
  Rutas creadas automáticamente. Repite por cada GSE/GSM en `data/processed/`.

---

### 4. `01.5.Filtrado.py`
- **Función**:  
  **Compilación de QC** y **barplots resumen** por GSE.
- **Incluye**:
  - Grids combinados por tipo de figura (violines, histogramas, scatter, top genes).
  - Barplots desde `qc_summary.csv`: `n_cells_final`, `mean_genes`, `mean_counts`, `mean_pct_mt`.
- **Salidas**:  
  `figures/01.filtered/<GSE>/figures_recopilation/*.png` + `<GSE>_qc_barplots_combined.png`.

---

### 5. `02.Normalization.py`
- **Función**:  
  **Normalización** (por total a 1e4 + `log1p`) y **selección de 2.000 HVGs** (Seurat).
- **Incluye**:
  - Procesado GSM-a-GSM desde `results/01.filtered/`.
  - Figuras HVG individuales y **grid combinado** por GSE.
- **Salidas**:
  - `results/02.normalized/<GSE>/<GSM>/normalized_data.h5ad`
  - `figures/02.normalized/<GSE>/<GSM>/highly_variable_genes_plot.png`
  - `figures/02.normalized/<GSE>/figures_recopilation/combined_hvg.png`

---

### 6. `03.Concatenation.py`
- **Función**:  
  **Concatenación multiestudio** (todas las muestras normalizadas).
- **Incluye**:
  - Añade metadatos `GSE`, `GSM`, `sample_id` en `.obs`.
  - Copia a `.raw` para usos posteriores (CellTypist).
- **Salidas**:
  - `results/03.concatenated/all_samples_concatenated.h5ad`
  - `results/03.concatenated/new_all_samples_concatenated_normalized_ready.h5ad`

---

### 7. `04.Integration_tutora.py`
- **Función**:  
  **Integración** con **Harmony** (por `GSM`) + **UMAP** pre/post + anotación de **condición clínica**.
- **Incluye**:
  - HVGs (Seurat v3 con `batch_key='GSM'`), `scale`, **PCA**, **Harmony** → `X_pca_harmony`.
  - Vecinos/UMAP pre/post; **UMAPs con marco negro**.
  - **Asignación de condición** por diccionario (`Sano`, `IPF`, `Ssc-ILD`, `COVID-19`).
  - **Submuestreo** controlado (hasta 50k células) para visualización global.
- **Salidas**:
  - `results/04.integration/all_samples_integrated_comparison.h5ad`
  - `results/04.integration/all_samples_umap_50k.h5ad`
  - `figures/04.integration/*umap*.png` (por GSE/condición/métricas)

---

### 8. `05.Clustering.py`
- **Función**:  
  **Clustering Leiden** y **selección de resolución** vía **estabilidad ARI**.
- **Incluye**:
  - 10 iteraciones sobre submuestreos (5.000 células) en resoluciones 0.3, 0.5 y 1.0.
  - Cálculo de **Adjusted Rand Index**; elección por **mediana ARI**.
  - Clustering final en el conjunto (espacio Harmony) y UMAP por clúster.
- **Salidas**:
  - `results/05.clustering/ari_scores.csv`
  - `results/05.clustering/all_samples_clustered_50k.h5ad`
  - `figures/05.clustering/ari_resolution_comparison.png`
  - `figures/05.clustering/umap_leiden_final_res_<X_X>.png`

---

### 9. `06.Final_analysis.py`
- **Función**:  
  **Anotación celular** (CellTypist por condición), **unificación de etiquetas**, **RGG** y **figuras de síntesis**.
- **Incluye**:
  - Modelos por condición:
    - Sano → `Human_Lung_Atlas.pkl`
    - COVID-19 → `Autopsy_COVID19_Lung.pkl`
    - IPF → `Human_IPF_Lung.pkl`
    - Ssc-ILD → `Human_PF_Lung.pkl`
  - Unificación de sinónimos (p. ej., “B”, “B cells” → “B cells”).
  - Filtro de tipos con `< 50` células → “Other”.
  - `rank_genes_groups` (Wilcoxon), **Heatmap top3 por Leiden**, **Dotplots top10**, **UMAPs** y **heatmap de proporciones**.
- **Salidas**:
  - `results/06.final_analysis/all_samples_annotated_celltype_unified.h5ad`
  - `figures/06.final_analysis/*` (UMAPs, heatmaps, dotplots, proporciones)

---

## Requisitos

- **Ramas y versiones**:
  - Python **≥ 3.9**
- **Paquetes**:
  - Core: `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`
  - Visualización: `matplotlib`, `seaborn`, `pillow`
  - Integración: `harmonypy`
  - ML: `scikit-learn`
  - Anotación: `celltypist`
- **Datos**:
  - Archivos crudos 10x (`matrix.h5`) **o** (`matrix.mtx.gz`, `genes/features.tsv.gz`, `barcodes.tsv.gz`) organizados en `data/processed/GSE*/GSM*/`.
  - Modelos de CellTypist accesibles (ver nombres arriba).

### Instalación rápida (conda + pip)
```bash
conda create -n tfm-scrna python=3.10 -y
conda activate tfm-scrna
pip install scanpy anndata harmonypy celltypist numpy pandas scipy matplotlib seaborn pillow scikit-learn
