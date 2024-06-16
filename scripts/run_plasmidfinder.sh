#!/bin/bash
#SBATCH --partition=longq7
#SBATCH -N 1
#SBATCH --mem-per-cpu=16000
#SBATCH --job-name=plasmidfinder
#SBATCH --output=plasmidfinder_%j.out 

# Load Singularity Module
module load singularity-3.8.7-gcc-9.4.0-xnxycrk

# Define the path for the Singularity Image File
plasmidfinder_sif="/mnt/beegfs/home/trane2023/containers/plasmidfinder.sif"

# Define the directory for the Plasmid Finder Database
plasmidfinder_db="/mnt/beegfs/home/trane2023/databases/plasmidfinder_db"

# Define the directory to the sample folders
input_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_assemblies"

# Define the directory that will contain the sample outputs
output_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_plasmid_data"
mkdir -p $output_dir

# Loop through each sample directory in the input directory
for sample_dir in "$input_dir"/*; do

    # Extract the name of the sample from the directory path
    sample_name=$(basename "$sample_dir")

    # Define the path to the genome assembly
    assembly_file=$(ls "$sample_dir/${sample_name}"*_genomic.fna)

    # Check if the assembly file exists
    if [ -f "$assembly_file" ]; then
        # Define the output directory for the sample
        sample_output_dir="$output_dir/$sample_name"
        
        # Create a directory for the sample in the output directory
        mkdir -p "$sample_output_dir"

        # Run PlasmidFinder using BLASTn
        singularity exec $plasmidfinder_sif plasmidfinder.py -i $assembly_file -o $sample_output_dir -p $plasmidfinder_db
    else
        echo "genomic.fna not found in $sample_dir"
    fi
done