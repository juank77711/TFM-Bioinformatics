#!/bin/bash

# Ir una carpeta atrás (asumimos que se ejecuta desde scripts/)
cd ..

# Definir rutas
RAW_DIR="data/raw"
PROCESSED_DIR="data/processed"

# Crear carpeta processed/ si no existe
mkdir -p "$PROCESSED_DIR"

# Entrar en raw
cd "$RAW_DIR" || { echo "❌ No se encontró $RAW_DIR"; exit 1; }

# Extraer los .tar si no están descomprimidos aún
for tar_file in GSE*_RAW.tar; do
    [[ -f "$tar_file" ]] || continue

    gse_name=$(echo "$tar_file" | cut -d'_' -f1)
    gse_dir="${gse_name}/"

    mkdir -p "$gse_dir"
    echo "📦 Descomprimiendo $tar_file en $gse_dir..."
    tar -xf "$tar_file" -C "$gse_dir"
done

# Reorganizar y renombrar archivos por GSM en data/processed/
for gse_dir in GSE*/; do
    gse_name="${gse_dir%/}"
    cd "$gse_dir" || continue

    for file in *; do
        [[ -f "$file" ]] || continue

        gsm_code=$(echo "$file" | grep -oE "GSM[0-9]+")
        [[ -z "$gsm_code" ]] && continue

        dest_dir="../../processed/$gse_name/$gsm_code"
        mkdir -p "$dest_dir"

        # 🧬 Si es .h5 filtrado, copiarlo como matrix.h5
        if [[ "$file" == *.h5 && "$file" == *filtered* ]]; then
            echo "🔹 Copiando H5 filtrado: $file → matrix.h5"
            cp "$file" "$dest_dir/matrix.h5"

        # 🧬 Si es barcodes
        elif [[ "$file" == *barcodes*.tsv.gz ]]; then
            echo "🔹 Copiando barcodes: $file → barcodes.tsv.gz"
            cp "$file" "$dest_dir/barcodes.tsv.gz"

        # 🧬 Si es genes o features
        elif [[ "$file" == *genes*.tsv.gz || "$file" == *features*.tsv.gz ]]; then
            echo "🔹 Copiando genes/features: $file → genes.tsv.gz"
            cp "$file" "$dest_dir/genes.tsv.gz"

        # 🧬 Si es la matriz
        elif [[ "$file" == *matrix*.mtx.gz ]]; then
            echo "🔹 Copiando matrix: $file → matrix.mtx.gz"
            cp "$file" "$dest_dir/matrix.mtx.gz"
        fi
    done

    cd ..
done

echo "✅ Archivos organizados y renombrados correctamente en data/processed/"

