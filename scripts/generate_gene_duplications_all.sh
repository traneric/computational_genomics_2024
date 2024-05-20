#!/bin/bash

# Navigate to the directory containing all the blast sample folders
cd /mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_blast_data
# Loop through each sample directory
for dir in */; do

    # Remove the trailing slash to get the directory name only
    sample_name="${dir%/}"

    # Define the path to the input results file
    input_file="${dir}results.txt"

    # Define the path for the output file
    output_file="${dir}${sample_name}_gene_duplications.txt"

    # Check if the results.txt file exists before attempting to process it
    if [[ -f "$input_file" ]]; then
        # Run awk command to filter and output to the sample-specific file
	# all-vs-all BLASTp, cutoff > 85%, alignment length between pairs > 85%, and eval < 10-10
        awk 'BEGIN{FS=OFS="\t"} $3 > 85 && (($8 - $7 + 1) / $8 >= 0.85) && $12 > 50 && $11 < 1e-10 && $1 != $2' "$input_file" > "$output_file"
        echo "Filtered data written to $output_file"
    else
        echo "Results file not found in $dir"
    fi
 done
