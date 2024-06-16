#!/bin/bash
#SBATCH --partition=longq7
#SBATCH -N 1
#SBATCH --mem-per-cpu=16000
#SBATCH --job-name=card_blast
#SBATCH --output=card_blast_%j.out 

# Activate your conda environment
source activate blast

# Define the input directory where the sample folders are located
input_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_assemblies"

# Define the output directory for BLAST results
output_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_amr_analysis_data"

# Define the CARD BLAST Database directory index
CARD_db="/mnt/beegfs/home/trane2023/computational_genomics_2024/card_database/card_protein_fasta_protein_homology_model/protein_db"

# Create the output directory if it doesn't already exist
mkdir -p "$output_dir"

# Loop through each sample directory in the input directory
for sample_dir in "$input_dir"/*; do
    # Extract the name of the sample from the directory path
    sample_name=$(basename "$sample_dir")

    # Define the path to the protein.faa file
    protein_file="$sample_dir/protein.faa"

    # Check if the protein.faa file exists
    if [ -f "$protein_file" ]; then
        # Define the output directory for the sample
        sample_output_dir="$output_dir/$sample_name"
        
        # Create a directory for the sample in the output directory
        mkdir -p "$sample_output_dir"

        # Define the output path for the BLAST results
        results_file="$sample_output_dir/results.txt"

        # Run BLASTp against the CARD BLAST DB, output the results to the sample's output directory
        blastp -query "$protein_file" -db "$CARD_db" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -out "$results_file"
    else
        echo "protein.faa not found in $sample_dir"
    fi
done

# Deactivate the conda environment
conda deactivate
