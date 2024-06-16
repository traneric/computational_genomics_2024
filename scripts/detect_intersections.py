#!/usr/bin/env python3

"""
The following script takes as input a GFF file, BLASTp output file in FMT 6, 
and a MobileElementFinder results file. It detects if the BLAST protein hits
reside inside or outside an MGE. It then calculates how many of ALL genes 
(duplicates and non-duplicates) reside inside or outside MGEs. The final output
is an MGE analysis report for the sample.
"""

import argparse
import re
from scipy.stats import binom_test

assembly_id = "GCA_024763995.1"
duplications_file = f"/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_gene_duplication_data/{assembly_id}_gene_duplications.txt"
gff_file = f"/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_assemblies/{assembly_id}/genomic.gff"
mge_file = f"/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_mobile_element_finder_data/{assembly_id}/mge.csv"
plasmid_file = f"/mnt/beegfs/home/trane2023/computational_genomics_2024/barnes_hospital_PRJNA824420/barnes_hospital_plasmid_data/{assembly_id}/plasmidfinder.csv"

def extract_duplicate_protein_coordinates(duplications_file, gff_file):
    coordinates_dict = {}
    subject_id_tracker = []

    with open(duplications_file, 'r') as dup_file, open(gff_file, 'r') as feature_file:
        for dup_line in dup_file:
            columns = dup_line.split('\t')
            subject_id = columns[1]

            # Since a subject_ID could appear more than once in the gene duplications file,
            # indicating that it is indeed a duplicate, there is no need to extract its 
            # coordinates again
            if subject_id in subject_id_tracker:
                continue
            subject_id_tracker.append(subject_id)

            # Reset the file pointer to the beginning of feature_file
            feature_file.seek(0)

            # Indentify the feature line in the GFF file that corresponds to this protein (subject-id)
            for feature_line in feature_file:
                if re.search(subject_id, feature_line):
                    feature_columns = feature_line.split('\t')
                    feature_contig_id = feature_columns[0]
                    start_coordinate = feature_columns[3]
                    end_coordinate = feature_columns[4]
                    duplicate_meta = [start_coordinate, end_coordinate, subject_id]
                 
                    if feature_contig_id not in coordinates_dict.keys(): 
                        coordinates_dict[feature_contig_id] = [duplicate_meta]
                    else:
                        coordinates_dict[feature_contig_id].append(duplicate_meta)

    coordinates_dict = combine_programmed_frameshift_coordinates(coordinates_dict)
    
    return coordinates_dict


def extract_all_protein_coordinates(gff_file):
    coordinates_dict = {}

    with open(gff_file, 'r') as gff_file:
        
        for line in gff_file:
            # Skip any header lines or blank lines
            if line.startswith('#') or not line.strip():
                continue
            
            # Split line into columns
            line = line.split('\t')
            feature = line[2]

            if feature == "CDS":
                sequence_region = line[0]
                sample_ID = line[8].split(';')[0]
                start_coordinate = line[3]
                end_coordinate = line[4]
                cds_meta = [start_coordinate, end_coordinate, sample_ID]

                if sequence_region not in coordinates_dict.keys():
                    coordinates_dict[sequence_region] = []
                coordinates_dict[sequence_region].append(cds_meta)


    # adjust the final coordinates to combine programmed frameshifts
    coordinates_dict = combine_programmed_frameshift_coordinates(coordinates_dict)

    return coordinates_dict


def combine_programmed_frameshift_coordinates(coordinates_dict):
    for contig_ID, tuples in coordinates_dict.items():
        subject_dict = {}
        for start, end, subject_id in tuples:
            if subject_id not in subject_dict:
                subject_dict[subject_id] = [start, end]
            else:
                subject_dict[subject_id][0] = min(subject_dict[subject_id][0], start)
                subject_dict[subject_id][1] = max(subject_dict[subject_id][1], end)
        coordinates_dict[contig_ID] = [(coords[0], coords[1], subject_id) for subject_id, coords in subject_dict.items()]
    return coordinates_dict


def extract_mge_coordinates(mge_file):
        mge_coordinate_dict = {}
        with open(mge_file, 'r') as mge_csv:

            # Skip the first six lines
            for _ in range(6):
                next(mge_csv)

            for mge_line in mge_csv:
                mge_cols = mge_line.split(',')
                mge_name = mge_cols[1]
                mge_prediction = mge_cols[3]
                mge_type = mge_cols[4]
                mge_contig_id = mge_cols[12].strip('"').split(" ")[0]
                mge_start_coordinate = mge_cols[14]
                mge_end_coordinate = mge_cols[15]
                mge_meta = (mge_start_coordinate, 
                            mge_end_coordinate,
                            mge_name,
                            mge_type)
                
                # Putative ComTns were not included in the MGE profile to avoid
                # introducing bias from false-positive or false-negative predictions
                # if (mge_prediction == "putative" and mge_type == "composite transposon"):
                #     continue

                if mge_contig_id not in mge_coordinate_dict.keys():
                    mge_coordinate_dict[mge_contig_id] = []
                mge_coordinate_dict[mge_contig_id].append(mge_meta)

        return mge_coordinate_dict


# def extract_plasmid_coordinates(plasmid_csv_file):
#     plasmid_coordinate_dict = {}

#     with open(plasmid_csv_file, 'r') as plasmid_csv:
#         for line in plasmid_csv:
#             print(line)


def detect_intersects(gene_coordinates, mge_coordinates):

    intersects = {} # Key :Contig_ID , Value: duplicate_data
    total_num_duplicates = sum(len(meta) for meta in gene_coordinates.values())
    num_duplicates_inside_mge = 0

    for contig_id, dup_meta_list in gene_coordinates.items():

        # No MGEs were found in this sequence_id
        if contig_id not in mge_coordinates.keys():
            continue
        
        # Initialize key to track duplicates that intersect in this sequence region
        if contig_id not in intersects.keys():
            intersects[contig_id] = []

        # Extract all mges that corresponds to this sequence_id
        sequence_mge_data = mge_coordinates[contig_id]

        # Iterate through each duplicate and see if it intersects with an MGE
        for dup_meta in dup_meta_list:
            dup_start = int(dup_meta[0]) - 31000
            dup_end = int(dup_meta[1]) + 31000

            for mge_meta in sequence_mge_data:
                mge_start = int(mge_meta[0])
                mge_end = int(mge_meta[1])

                # Check if the ranges intersect
                if (dup_start <= mge_end and dup_end >= mge_start):
                    intersects[contig_id].append(dup_meta)
                    num_duplicates_inside_mge += 1
                    break

    num_duplicates_outside_mge = total_num_duplicates - num_duplicates_inside_mge

    intersection_analysis_data = [total_num_duplicates, num_duplicates_inside_mge, num_duplicates_outside_mge]

    return intersection_analysis_data


def generate_final_report(assembly_id, dup_mge_data, all_mge_data):
    
    total_num_genes = all_mge_data[0]
    total_num_genes_inside_mge = all_mge_data[1]
    total_num_genes_outside_mge = all_mge_data[2]
    total_genes_intersection_percentage = round((total_num_genes_inside_mge / total_num_genes) * 100, 2)

    num_duplicate_genes = dup_mge_data[0]
    num_duplicate_genes_inside_mge = dup_mge_data[1]
    num_duplicate_genes_outside_mge = dup_mge_data[2]
    duplicate_genes_intersection_percentage = round((num_duplicate_genes_inside_mge / num_duplicate_genes) * 100, 2)
    
    x = num_duplicate_genes_inside_mge # num successes
    n = total_num_genes_inside_mge # total number of trials
    p = num_duplicate_genes_outside_mge / total_num_genes_outside_mge # the probability of success under the null hypothesis
    p_value = binom_test(x, n, p, alternative='two-sided')
    p_value = format(p_value, ".2e")

    # Print results
    print(assembly_id,
          total_num_genes, 
          total_num_genes_inside_mge, 
          total_num_genes_outside_mge, 
          total_genes_intersection_percentage, 
          num_duplicate_genes, 
          num_duplicate_genes_inside_mge, 
          num_duplicate_genes_outside_mge,
          duplicate_genes_intersection_percentage,
          p_value, sep=",")


def main(assembly_id, duplications_file, gff_file, mge_file):
    mge_coordinates = extract_mge_coordinates(mge_file)
    # plasmid_coordinates = extract_plasmid_coordinates(plasmid_file)

    dup_protein_coordinates = extract_duplicate_protein_coordinates(duplications_file, gff_file)
    dup_mge_analysis_data = detect_intersects(dup_protein_coordinates, mge_coordinates)

    all_protein_coordinates = extract_all_protein_coordinates(gff_file)
    all_mge_analysis_data = detect_intersects(all_protein_coordinates, mge_coordinates)

    generate_final_report(assembly_id, dup_mge_analysis_data, all_mge_analysis_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process gene files")
    parser.add_argument("--assembly_id", type=str, help="Assembly ID")
    parser.add_argument("--duplications_file", type=str, help="Path to the gene duplications file")
    parser.add_argument("--gff_file", type=str, help="Path to the GFF file")
    parser.add_argument("--mge_file", type=str, help="Path to the MGE file")

    args = parser.parse_args()

    main(args.assembly_id, args.duplications_file, args.gff_file, args.mge_file)
# main(assembly_id, duplications_file, gff_file, mge_file, plasmid_file)
