#!/bin/bash

conda activate python

# Define the path to the python script
python_script="/mnt/beegfs/home/trane2023/computational_genomics_2024/scripts/generate_plasmid_csv.py"

# Define the directory to the PlasmidFinder sample folders
input_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_plasmid_data"

# Loop through each sample directory in the input directory
for sample_dir in "$input_dir"/*; do

    # Define the path to the PlasmidFinder JSON File
    json_file="$sample_dir/data.json"

    # Check if the JSON file exists
    if [ -f "$json_file" ]; then
        # Define the output CSV file
        echo "processing $json_file"
        csv_output_file="$sample_dir/plasmidfinder.csv"
        
        python $python_script --json_path $json_file >> $csv_output_file
    else
        echo "data.json not found in $sample_dir"
    fi
done