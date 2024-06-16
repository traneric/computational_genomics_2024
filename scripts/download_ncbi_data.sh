#!/bin/bash
# Directory where you want to store all the downloaded data
PARENT_DIR="/path/to/downloaded_data"

# File containing the list of accession numbers
ACCESSION_FILE="accessions.txt"

# Create parent directory if it doesn't exist
mkdir -p "$PARENT_DIR"

# Read each accession number and process it
while IFS= read -r ACCESSION; do

    # Create a unique directory for each accession
    ACCESSION_DIR="${PARENT_DIR}/${ACCESSION}"
    mkdir -p "$ACCESSION_DIR"

    # Change to the directory for the accession
    cd "$ACCESSION_DIR"

    # Download data for the accession
    datasets download genome accession "$ACCESSION" --include gff3,rna,cds,protein,genome,seq-report --filename "${ACCESSION}.zip"
   
    # Unzip the downloaded data
    unzip "${ACCESSION}.zip"

    # Optionally, remove the zip file to save space
    rm "${ACCESSION}.zip"

    # Change back to the parent directory
    cd "$PARENT_DIR" done < "$ACCESSION_FILE" echo "All downloads completed."
