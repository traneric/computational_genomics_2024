#!/bin/bash
#SBATCH --partition=mediumq7
#SBATCH -N 1
#SBATCH --mem-per-cpu=16000
#SBATCH --job-name=mobile_element_finder
#SBATCH --output=mobile_element_finder%j.out 

# Activate your conda environment
source activate MobileElementFinder

# Define the input directory where the sample folders are located
PARENT_DIR="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_assemblies"

# Output directory for MobileElementFinder results
OUTPUT_DIR="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_mobile_element_finder_data"

# Loop through each sample directory
for SAMPLE_DIR in "$PARENT_DIR"/*; do
  # Check if it is a directory
  if [ -d "$SAMPLE_DIR" ]; then
    # Extract sample name
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")

    # Find the genomic.fna file
    GENOMIC_FILE=$(find "$SAMPLE_DIR" -type f -name "${SAMPLE_NAME}*_genomic.fna")
    
    # Check if the genomic.fna file exists
    if [ -f "$GENOMIC_FILE" ]; then
      SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_NAME"
      mkdir -p "$SAMPLE_OUTPUT_DIR"
      
      # Run MGE Finder with the genomic.fna file
      echo "Running mefinder on $GENOMIC_FILE for sample $SAMPLE_NAME"
      mefinder find --contig $GENOMIC_FILE --gff mge
      mv mge* "$SAMPLE_OUTPUT_DIR/"
    else
      echo "No genomic.fna file found in $SAMPLE_DIR"
    fi
  fi 
done

