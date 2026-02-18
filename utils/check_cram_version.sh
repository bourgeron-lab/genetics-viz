#!/bin/bash
# check_cram_version.sh - Check CRAM file version for IGV.js compatibility
#
# IGV.js uses the cram-js library which supports CRAM v2.0, v2.1, and v3.0.
# CRAM v3.1 introduces codecs (adaptive arithmetic, name tokenizer, fqzcomp)
# that are not supported. Files needing recompression can be fixed with
# recompress_cram.sh.
#
# Exit codes:
#   0 - CRAM file is compatible with IGV.js
#   1 - Error (bad arguments, file not found, invalid file)
#   2 - CRAM file needs recompression for IGV.js compatibility

set -e
set -u

# --- Argument validation ---
if [ $# -eq 0 ]; then
    echo "Error: Please provide a barcode as parameter" >&2
    echo "Usage: $0 <BARCODE>" >&2
    exit 1
fi

BARCODE=$1

# Validate barcode format (digits and capital letters only)
if ! [[ "$BARCODE" =~ ^[A-Z0-9]+$ ]]; then
    echo "Error: Barcode must contain only digits and capital letters" >&2
    exit 1
fi

# --- Path setup ---
BASE_PATH="/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/samples/${BARCODE}/sequences"
CRAM_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram"
CRAI_FILE="${BASE_PATH}/${BARCODE}.GRCh38_GIABv3.cram.crai"

# --- File checks ---
if [ ! -f "$CRAM_FILE" ]; then
    echo "Error: CRAM file not found: $CRAM_FILE" >&2
    exit 1
fi

if [ ! -s "$CRAM_FILE" ]; then
    echo "Error: CRAM file is empty: $CRAM_FILE" >&2
    exit 1
fi

# --- Read CRAM version from magic bytes ---
# CRAM format: bytes 1-4 = "CRAM" (magic), byte 5 = major version, byte 6 = minor version
HEADER_HEX=$(xxd -l 6 -p "$CRAM_FILE")

# Validate magic number (CRAM = 43 52 41 4d in hex)
MAGIC="${HEADER_HEX:0:8}"
if [ "$MAGIC" != "4352414d" ]; then
    echo "Error: Not a valid CRAM file (bad magic number): $CRAM_FILE" >&2
    exit 1
fi

# Extract version (hex -> decimal)
MAJOR=$((16#${HEADER_HEX:8:2}))
MINOR=$((16#${HEADER_HEX:10:2}))
VERSION="${MAJOR}.${MINOR}"

# --- Gather additional info ---
FILE_SIZE=$(du -h "$CRAM_FILE" | cut -f1)
CRAI_STATUS="Not found"
if [ -f "$CRAI_FILE" ]; then
    CRAI_STATUS="Found"
fi

# --- Output report ---
echo "=== CRAM Version Check ==="
echo "Barcode:       $BARCODE"
echo "CRAM file:     $CRAM_FILE"
echo "File size:     $FILE_SIZE"
echo "CRAI index:    $CRAI_STATUS"
echo "CRAM version:  $VERSION"
echo ""

# --- Compatibility decision ---
# Whitelist of versions known to work with IGV.js (cram-js)
COMPATIBLE_VERSIONS=("2.0" "2.1" "3.0")

IS_COMPATIBLE=false
for v in "${COMPATIBLE_VERSIONS[@]}"; do
    if [ "$VERSION" = "$v" ]; then
        IS_COMPATIBLE=true
        break
    fi
done

if $IS_COMPATIBLE; then
    echo "RESULT: COMPATIBLE"
    echo "  CRAM version $VERSION is supported by IGV.js."
    if [ "$CRAI_STATUS" = "Not found" ]; then
        echo ""
        echo "  WARNING: CRAI index file is missing. IGV.js requires an index."
        echo "  Run: samtools index $CRAM_FILE"
    fi
    exit 0
else
    echo "RESULT: NEEDS RECOMPRESSION"
    echo "  CRAM version $VERSION uses codecs not supported by IGV.js (cram-js)."
    echo "  Run: ./recompress_cram.sh $BARCODE"
    exit 2
fi
