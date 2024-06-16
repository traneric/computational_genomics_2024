#!/bin/bash

# Define the parent directory containing subdirectories with results.txt files
parent_directory="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_amr_analysis_data"

# Function to process each results.txt file
process_results() {
    local input_file="$1/results.txt"
    local output_file="$1/top_amr_hits.txt"
    
# Check if the results.txt file exists
    if [[ -f $input_file ]]; then
        # Step 1: Filter the BLASTp results based on the given criteria
        awk '$11 <= 1e-10 && $3 >= 90 && $4/($4+$6) >= 0.5 && $12 >= 50' $input_file > filtered_tmp.txt
        
	# Step 2: Sort by query ID (column 1) and bit score (column 12) in descending order
        sort -k1,1 -k12,12nr filtered_tmp.txt > sorted_filtered_tmp.txt
        
	# Step 3: Use awk to keep only the best hit for each query
        awk '!seen[$1]++' sorted_filtered_tmp.txt > $output_file
        
	# Clean up temporary files
        rm filtered_tmp.txt sorted_filtered_tmp.txt
        
	echo "Filtering and selection of best hits completed for $input_file. Results saved to $output_file."
    else
        echo "No results.txt file found in $1"
    fi
}

# Loop through each subdirectory in the parent directory
for dir in "$parent_directory"/*/; do
    # Call the function to process results.txt in each subdirectory
    process_results "$dir" 
done
