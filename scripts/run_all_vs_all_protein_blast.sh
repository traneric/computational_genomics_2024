#!/bin/bash
#SBATCH --partition=mediumq7
#SBATCH -N 1
#SBATCH --mem-per-cpu=16000
#SBATCH --job-name=blast_job
#SBATCH --output=blast_job_%j.out 

# Activate your conda environment
source activate blast

# Define the input directory where the sample folders are located
input_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_assemblies"

# Define the output directory for BLAST results
output_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_blast_data"

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

        # Run makeblastdb, direct the database files to the sample's output directory
        makeblastdb -in "$protein_file" -dbtype prot -out "$sample_output_dir/protein_db"

        # Run BLASTp, output the results to the sample's output directory
        blastp -query "$protein_file" -db "$sample_output_dir/protein_db" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -out "$results_file"
    else
        echo "protein.faa not found in $sample_dir"
    fi
done

# Deactivate the conda environment
conda deactivate
