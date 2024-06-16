#!/bin/bash

# Directories
# GFF_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_assemblies"
# duplications_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_gene_duplication_data"
# MGE_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_mobile_element_finder_data"

GFF_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_assemblies"
duplications_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_gene_duplication_data"
MGE_dir="/mnt/beegfs/home/trane2023/computational_genomics_2024/carbapenem_res_PRJEB30134/PRJEB30134_mobile_element_finder_data"

# Create CSV Header
echo "assembly_id,total_num_genes,total_num_genes_inside_mges,total_num_genes_outside_mges,percent_total_genes_in_mge,num_duplicate_genes,num_duplicate_genes_inside_mge,num_duplicate_genes_outside_mge,percent_duplicate_genes_in_mge,p_val" >> barnes_duplicate_mge_report_new.csv

# Loop over each assembly directory in GFF_dir
for assembly_dir in "$GFF_dir"/*/; do
    # Extract the assembly ID from the directory path
    assembly_id=$(basename "$assembly_dir")

    # Define the path to the genomic.gff file
    gff_file="${assembly_dir}genomic.gff"

    # Define the path to the gene duplications file
    duplications_file="${duplications_dir}/${assembly_id}_gene_duplications.txt"

    # Define the path to the MGE file
    mge_file="${MGE_dir}/${assembly_id}/mge.csv"

    # Check if the necessary files exist
    if [ -f "$gff_file" ] && [ -f "$duplications_file" ] && [ -f "$mge_file" ]; then
        echo "Processing $assembly_id"
        echo "GFF File: $gff_file"
        echo "Gene Duplications File: $duplications_file"
        echo "MGE File: $mge_file"
        
        # Call your Python script here with paths to GFF, gene duplications, and MGE files
        python detect_intersections.py --assembly_id $assembly_id --gff_file $gff_file --duplications_file $duplications_file --mge_file $mge_file >> carbapenem_duplicate_mge_report_ext.csv
    else
        echo "Required files not found for $assembly_id"
        [ ! -f "$gff_file" ] && echo "Missing: $gff_file"
        [ ! -f "$duplications_file" ] && echo "Missing: $duplications_file"
        [ ! -f "$mge_file" ] && echo "Missing: $mge_file"
    fi
done
