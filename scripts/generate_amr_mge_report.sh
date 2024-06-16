#!/bin/bash

# Directories
GFF_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_assemblies"
amr_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_amr_analysis_data/CARD_blast_output"
MGE_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_mobile_element_finder_data"

# GFF_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_assemblies"
# amr_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_amr_analysis_data"
# MGE_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_mobile_element_finder_data"

# Loop over each assembly directory in GFF_dir
for assembly_dir in "$GFF_dir"/*/; do
    # Extract the assembly ID from the directory path
    assembly_id=$(basename "$assembly_dir")

    # Define the path to the genomic.gff file
    gff_file="${assembly_dir}genomic.gff"

    # Define the path to the amr top hits file
    amr_file="${amr_dir}/${assembly_id}/top_amr_hits.txt"

    # Define the path to the MGE file
    mge_file="${MGE_dir}/${assembly_id}/mge.csv"

    # Check if the necessary files exist
    if [ -f "$gff_file" ] && [ -f "$amr_file" ] && [ -f "$mge_file" ]; then
        echo "Processing $assembly_id"
        echo "GFF File: $gff_file"
        echo "AMR Top Hits File: $amr_file"
        echo "MGE File: $mge_file"
        
        # Call your Python script here with paths to GFF, gene duplications, and MGE files
        python detect_amr_mge_intersections.py --assembly_id $assembly_id --gff_file $gff_file --amr_file $amr_file --mge_file $mge_file >> carbapenem_amr_mge_report_ext.csv
    else
        echo "Required files not found for $assembly_id"
        [ ! -f "$gff_file" ] && echo "Missing: $gff_file"
        [ ! -f "$amr_file" ] && echo "Missing: $amr_file"
        [ ! -f "$mge_file" ] && echo "Missing: $mge_file"
    fi
done
