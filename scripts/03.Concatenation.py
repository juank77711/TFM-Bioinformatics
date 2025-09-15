import os
from pathlib import Path
import scanpy as sc
import anndata as ad

# Rutas base
base_dir = Path(__file__).resolve().parent.parent
normalized_base = base_dir / "results" / "02.normalized"
concat_dir = base_dir / "results" / "03.concatenated"
concat_dir.mkdir(parents=True, exist_ok=True)

# Ruta del archivo final concatenado
output_path = concat_dir / "all_samples_concatenated.h5ad"

adatas = []

# Recorremos todas las carpetas GSE dentro de 02.normalized
for gse_dir in normalized_base.iterdir():
    if not gse_dir.is_dir() or not gse_dir.name.startswith("GSE"):
        continue

    for gsm_dir in gse_dir.iterdir():
        if not gsm_dir.is_dir() or not gsm_dir.name.startswith("GSM"):
            continue

        h5ad_path = gsm_dir / "normalized_data.h5ad"
        if h5ad_path.exists():
            print(f"🔄 Cargando {gse_dir.name}/{gsm_dir.name}...")
            adata = sc.read_h5ad(h5ad_path)

            # Añadir metadatos
            adata.obs["GSE"] = gse_dir.name
            adata.obs["GSM"] = gsm_dir.name
            adatas.append(adata)

# Concatenar si hay datos
if len(adatas) == 0:
    print("⚠️ No se encontraron datos para concatenar.")
else:
    print("🧬 Concatenando todos los objetos AnnData...")
    adata_concat = ad.concat(adatas, join='outer', label='sample_id', keys=[a.obs['GSM'][0] for a in adatas])

    # Guardar el archivo concatenado original
    adata_concat.write(output_path)
    print(f"✅ Archivo concatenado guardado en: {output_path}")

    # Crear y guardar el preparado para CellTypist
    adata_concat.raw = adata_concat.copy()
    ready_path = concat_dir / "new_all_samples_concatenated_normalized_ready.h5ad"
    adata_concat.write(ready_path)
    print(f"📦 `.raw` creado y guardado en: {ready_path}")

