#!/bin/bash
# Migration script for snvs.tsv validation file
# Migrates from version 0.1.1 format to 0.2.0 format
#
# Old format (0.1.1):
#   FID	Variant	Sample	User	Inheritance	Validation	Timestamp
#
# New format (0.2.0):
#   FID	Variant	Sample	User	Inheritance	Validation	Comment	Ignore	Timestamp
#
# Changes:
#   - Added "Comment" column (empty string default)
#   - Added "Ignore" column (0 default = not ignored)
#
# Usage:
#   ./snvs_validations_migration_0.1.1_to_0.2.0.sh /path/to/snvs.tsv

set -e

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 /path/to/snvs.tsv"
    echo ""
    echo "Migrates snvs.tsv from version 0.1.1 to 0.2.0 format."
    echo "A backup will be created as snvs.tsv.backup"
    exit 1
fi

INPUT_FILE="$1"

# Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File not found: $INPUT_FILE"
    exit 1
fi

# Check if file is readable
if [ ! -r "$INPUT_FILE" ]; then
    echo "Error: Cannot read file: $INPUT_FILE"
    exit 1
fi

# Read the header
HEADER=$(head -n 1 "$INPUT_FILE")

# Check if already migrated (has Comment and Ignore columns)
if echo "$HEADER" | grep -q $'\tComment\t' && echo "$HEADER" | grep -q $'\tIgnore\t'; then
    echo "File appears to already be in 0.2.0 format (has Comment and Ignore columns)."
    echo "No migration needed."
    exit 0
fi

# Check if it's the expected old format
EXPECTED_OLD_HEADER="FID	Variant	Sample	User	Inheritance	Validation	Timestamp"
if [ "$HEADER" != "$EXPECTED_OLD_HEADER" ]; then
    echo "Error: Unexpected header format."
    echo "Expected: $EXPECTED_OLD_HEADER"
    echo "Found:    $HEADER"
    echo ""
    echo "This script only migrates from 0.1.1 to 0.2.0 format."
    exit 1
fi

# Create backup
BACKUP_FILE="${INPUT_FILE}.backup"
if [ -f "$BACKUP_FILE" ]; then
    echo "Backup file already exists: $BACKUP_FILE"
    read -p "Overwrite backup? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 1
    fi
fi

cp "$INPUT_FILE" "$BACKUP_FILE"
echo "Created backup: $BACKUP_FILE"

# Create new header
NEW_HEADER="FID	Variant	Sample	User	Inheritance	Validation	Comment	Ignore	Timestamp"

# Create temporary file for the migration
TEMP_FILE=$(mktemp)

# Write new header
echo "$NEW_HEADER" > "$TEMP_FILE"

# Process data rows (skip header)
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r FID Variant Sample User Inheritance Validation Timestamp; do
    # Add empty Comment and 0 for Ignore between Validation and Timestamp
    echo -e "${FID}\t${Variant}\t${Sample}\t${User}\t${Inheritance}\t${Validation}\t\t0\t${Timestamp}" >> "$TEMP_FILE"
done

# Replace original file
mv "$TEMP_FILE" "$INPUT_FILE"

echo "Migration complete!"
echo ""
echo "Old format: FID|Variant|Sample|User|Inheritance|Validation|Timestamp"
echo "New format: FID|Variant|Sample|User|Inheritance|Validation|Comment|Ignore|Timestamp"
echo ""
echo "All existing rows now have:"
echo "  - Comment: (empty)"
echo "  - Ignore: 0 (not ignored)"
