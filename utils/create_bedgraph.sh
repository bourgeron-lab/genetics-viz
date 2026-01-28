#!/bin/bash

# Script to generate bedgraph files from CRAM using mosdepth
# Usage: ./create_bedgraph.sh <BARCODE> [BIN_SIZE]

# Check if barcode parameter is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide a barcode as parameter"
    echo "Usage: $0 <BARCODE> [BIN_SIZE]"
    echo "  BARCODE  - Sample barcode (required)"
    echo "  BIN_SIZE - Bin size in bp (optional, default: 1000)"
    exit 1
fi

BARCODE=$1
BIN_SIZE=${2:-1000}

# Validate barcode format (digits and capital letters only)
if ! [[ "$BARCODE" =~ ^[A-Z0-9]+$ ]]; then
    echo "Error: Barcode must contain only digits and capital letters"
    exit 1
fi

# Validate bin size is a positive integer
if ! [[ "$BIN_SIZE" =~ ^[0-9]+$ ]] || [ "$BIN_SIZE" -eq 0 ]; then
    echo "Error: Bin size must be a positive integer"
    exit 1
fi

# Define base path for cleaner code
BASE_PATH="/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/samples/${BARCODE}/sequences"
CRAM_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram"

# Check if source CRAM file exists
if [ ! -f "$CRAM_FILE" ]; then
    echo "Error: CRAM file not found: $CRAM_FILE"
    exit 1
fi

# Create SLURM script content - Note the quoted 'EOF' to prevent variable expansion
cat <<'EOF' | sed "s|{{BARCODE}}|${BARCODE}|g; s|{{BIN_SIZE}}|${BIN_SIZE}|g" | sbatch
#!/bin/bash
#SBATCH --job-name=create_bedgraph_{{BARCODE}}
#SBATCH --partition=ghfc
#SBATCH --qos=ghfc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -e  # Exit on error
set -u  # Exit on undefined variable

BARCODE="{{BARCODE}}"
BIN_SIZE="{{BIN_SIZE}}"
BASE_PATH="/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/samples/${BARCODE}/sequences"
CRAM_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram"
REF_FASTA="/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta"
MOSDEPTH="/pasteur/helix/projects/ghfc_wgs/tools/mosdepth-0.3.12/mosdepth"
OUTPUT_PREFIX="${BASE_PATH}/${BARCODE}.by${BIN_SIZE}"

echo "Starting bedgraph generation for barcode: ${BARCODE}"
echo "Bin size: ${BIN_SIZE}"
echo "Job started at: $(date)"

# Change to output directory
cd "${BASE_PATH}"

# Limit virtual memory to prevent silent crashes from graceful degradation
ulimit -v $((8 * 1024 * 1024))  # 8GB in KB

# Run mosdepth
# -t: number of threads
# -b: bin size
# -n: don't output per-base depth (faster)
# -f: reference fasta
# -x: don't look at internal cigar operations or determine fragment coordinates based on mate position
# -Q: mapping quality threshold (40)
# -m: median mode (for CNV calling)
echo "Running mosdepth..."
${MOSDEPTH} -t 2 -b ${BIN_SIZE} -n -f ${REF_FASTA} -x -Q 40 -m ${OUTPUT_PREFIX} ${CRAM_FILE}

echo "Mosdepth completed at $(date)"

# Wait for filesystem to sync
sleep 5

# Rename the output file to match expected naming
echo "Renaming output file..."
mv ${OUTPUT_PREFIX}.regions.bed.gz ${OUTPUT_PREFIX}.bedgraph.gz

# Also rename the index file if it exists
if [ -f "${OUTPUT_PREFIX}.regions.bed.gz.csi" ]; then
    mv ${OUTPUT_PREFIX}.regions.bed.gz.csi ${OUTPUT_PREFIX}.bedgraph.gz.csi
fi

# Create tabix index for the bedgraph
echo "Creating tabix index..."
module load htslib/1.21 2>/dev/null || true
tabix -p bed ${OUTPUT_PREFIX}.bedgraph.gz

# Clean up mosdepth auxiliary files
echo "Cleaning up auxiliary files..."
rm -f ${OUTPUT_PREFIX}.mosdepth.global.dist.txt
rm -f ${OUTPUT_PREFIX}.mosdepth.region.dist.txt
rm -f ${OUTPUT_PREFIX}.mosdepth.summary.txt

echo "Bedgraph generation completed successfully for barcode: ${BARCODE}"
echo "Output files:"
echo "  - ${OUTPUT_PREFIX}.bedgraph.gz"
echo "  - ${OUTPUT_PREFIX}.bedgraph.gz.tbi"
echo "Job finished at: $(date)"
EOF

echo "Job submitted for barcode: ${BARCODE} with bin size: ${BIN_SIZE}"
