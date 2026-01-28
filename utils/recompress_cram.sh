#!/bin/bash

# Check if barcode parameter is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide a barcode as parameter"
    echo "Usage: $0 <BARCODE>"
    exit 1
fi

BARCODE=$1

# Validate barcode format (digits and capital letters only)
if ! [[ "$BARCODE" =~ ^[A-Z0-9]+$ ]]; then
    echo "Error: Barcode must contain only digits and capital letters"
    exit 1
fi

# Define base path for cleaner code
BASE_PATH="/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/samples/${BARCODE}/sequences"
CRAM_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram"
CRAI_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram.crai"

# Check if source files exist
if [ ! -f "$CRAM_FILE" ]; then
    echo "Error: CRAM file not found: $CRAM_FILE"
    exit 1
fi

# Create SLURM script content - Note the quoted 'EOF' to prevent variable expansion
cat <<'EOF' | sed "s|{{BARCODE}}|${BARCODE}|g" | sbatch
#!/bin/bash
#SBATCH --job-name=cram_recompression
#SBATCH --partition=ghfc
#SBATCH --qos=ghfc
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=cram_recompression_{{BARCODE}}_%j.out
#SBATCH --error=cram_recompression_{{BARCODE}}_%j.err

set -e  # Exit on error
set -u  # Exit on undefined variable

BARCODE="{{BARCODE}}"
BASE_PATH="/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/samples/${BARCODE}/sequences"
CRAM_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram"
CRAI_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram.crai"
PROBLEM_CRAM="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.problem.cram"
PROBLEM_CRAI="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.problem.cram.crai"
REF_FASTA="/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta"

echo "Starting CRAM recompression for barcode: ${BARCODE}"
echo "Job started at: $(date)"

# 1. Load samtools module
module load samtools/1.21

# 2. Move original files to .problem versions
echo "Moving original CRAM and CRAI files to .problem versions..."
mv "${CRAM_FILE}" "${PROBLEM_CRAM}"
if [ -f "${CRAI_FILE}" ]; then
    mv "${CRAI_FILE}" "${PROBLEM_CRAI}"
fi

# 3. Recompress CRAM file
echo "Recompressing CRAM file..."
samtools view -C -T "${REF_FASTA}" -O cram,level=1 "${PROBLEM_CRAM}" -o "${CRAM_FILE}"

# 4. Index the new CRAM file
echo "Indexing new CRAM file..."
samtools index "${CRAM_FILE}"

# 5. Remove .problem files
echo "Removing .problem files..."
rm -f "${PROBLEM_CRAM}"
rm -f "${PROBLEM_CRAI}"

echo "CRAM recompression completed successfully for barcode: ${BARCODE}"
echo "Job finished at: $(date)"
EOF

echo "Job submitted for barcode: ${BARCODE}"