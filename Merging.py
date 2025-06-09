import gffpandas.gffpandas as gffpd
import pandas as pd

import numpy as np
import os

# raw path for PC
raw_path = 'D:/Studies/MSc/First year/2 - Second semester/Research methods in Bioinformatics/Final_Project/GTEx_Data/'

# # raw path for laptop
# raw_path = 'C:/Users/eliro/Desktop/Studies/MSc/First year/2 - Second semester/Research methods in Bioinformatics/Final_Project/GTEx_Data/'

def separate_age_groups():
    # @title Making seperate files for each age group
    # tpm_file_path = data_raw_path + 'Transcripts_TPM_GTEx_Analysis_v10_RSEMv1.3.3.txt.gz'
    tpm_file_path = raw_path + 'RawData_v10/Transcripts_TPM_GTEx_Analysis_v10_RSEMv1.3.3.txt.gz'

    # Loading the sample_attributes_phenotypes file
    subj_path = raw_path + 'Output Files/Sample_Attributes_Phenotypes.txt'
    subj_att_pheno = pd.read_csv(subj_path, delimiter='\t')

    # age_list = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79']
    age_list = ['50-59', '60-69', '70-79']
    ages = dict()

    for agerange in age_list:
      # Extracting the IDs of the samples for the current age group range (E.g, 20-29)
      ages[agerange] = subj_att_pheno[(subj_att_pheno['AGE'] == agerange)].SAMPID.tolist()

      # The list that contains the names of the columns to extract for each age group
      list_of_cols = ['transcript_id','gene_id'] + ages[agerange]

      # Reads from the data a specific age group
      specific_age_group = pd.read_csv(tpm_file_path, delimiter='\t', usecols=list_of_cols)

      # Saving each age group in a different file
      file_name = raw_path + 'Output Files/' + f"Age Groups/{agerange}_Transcripts_TPM.csv"
      # file_name = data_raw_path + f"Age Groups/{agerange}_Transcripts_TPM.csv"
      specific_age_group.to_csv(file_name, sep='\t', index=False)
      print()


def loading_annotation():
    # @title Loading the gene annotation files, and basic filtering

    """
    In this block, I'm loading the gene annotation file.
    Right now, I'm using only the 3 annotations - start, end and transcript_id:
    - start: The start position of the transcript.
    - end: The end position of the transcript.
    - transcript_id: The ID of the transcript.

    *Note - there are many more annotations to use, but I'm using only these 3 for now.
    *Note - the origin of the annotation file is GENCODE - https://www.gencodegenes.org/human/.
    """

    raw_data_path = raw_path + 'gencode.v48.annotation.gff3'
    annotation = gffpd.read_gff3(raw_data_path)

    # # Deleting all the numbers after the decimal number
    acc_n = annotation.filter_feature_of_type(['gene']).attributes_to_columns().gene_id.str.split('.', expand=True)[0]

    # Choosing the columns I want
    gene_annotation = annotation.filter_feature_of_type(['transcript']).attributes_to_columns()[
        ['transcript_id', 'start', 'end', 'strand']].join(acc_n)

    # Cleaning the data - in 'transcript_id' column, deleting the decimal number and all the numbers after it
    gene_annotation['transcript_id'] = gene_annotation['transcript_id'].str.replace(r'\.[A-Za-z0-9]+$', '', regex=True)

    # Remove the last column of df - because it's just filled with NaNs
    gene_annotation.drop(gene_annotation.columns[-1], axis=1, inplace=True)

    # Adding 2 new columns which will be the start/end site according to the strand (+ or -)
    gene_annotation['true_start'] = np.where(  # Adding true_start column
        gene_annotation['strand'] == '+',
        gene_annotation['start'],
        gene_annotation['end']
    )

    # Dropping the columns strand, end and start
    gene_annotation.drop(columns=['strand', 'end', 'start'], axis=1, inplace=True)

    return gene_annotation


def merging():
    gene_annotation = loading_annotation()

    age_dir_path = raw_path + 'Output Files/Age Groups/'

    for age_group_file_name in os.listdir(age_dir_path):
        # Processing only files that ends with .csv
        if age_group_file_name.endswith('.csv'):
            file_path = os.path.join(age_dir_path, age_group_file_name)

            age_group_df = pd.read_csv(file_path, delimiter='\t')

            # Cleaning - deleting numbers after the decimal point
            for col in age_group_df.columns[:2]:
                age_group_df[col] = age_group_df[col].str.replace(r'\.[A-Za-z0-9]+$', '', regex=True)

            # Merging the two dfs based on the 'transcript_id' column
            merged_df = pd.merge(gene_annotation, age_group_df, on='transcript_id', how='inner')

            merged_df_file_path = age_dir_path + 'Age Groups With true_start/' + f"{age_group_file_name}"
            merged_df.to_csv(merged_df_file_path, sep='\t', index=False)
            print(f"Done saving the file {merged_df_file_path}")

            # After merging, dropping the 'transcript_id' column - no more use
            merged_df.drop(columns=['transcript_id'], axis=1, inplace=True)


            # Creating Two Tables according to TSS - looking at the pair values of 'true_start' and 'gene_id'
            # Mark duplicates based on both columns
            duplicate_mask = merged_df.duplicated(subset=['true_start', 'gene_id'], keep=False)

            # First table - NON-unique (true_start, gene_id) pairs
            df_duplicates = merged_df[duplicate_mask].copy()
            # Grouping rows based on (true_start, gene_id) pairs
            df_duplicates_combined = df_duplicates.groupby(['true_start', 'gene_id'], as_index=False).sum()
            df_duplicates_combined.sort_values(by='gene_id', inplace=True)

            # Second table - unique (true_start, gene_id) pairs
            df_unique = merged_df[~duplicate_mask].copy()


            # Saving the NON uniques after summing
            non_unique_tss_file_path = age_dir_path + 'NON_UniqueTSS/' + f"NON_UniqueTSS_{age_group_file_name}"
            df_duplicates_combined.to_csv(non_unique_tss_file_path, sep='\t', index=False)
            print(f"Done saving the file {non_unique_tss_file_path}")

            # Saving the uniques
            unique_tss_file_path = age_dir_path + 'UniqueTSS/' + f"UniqueTSS_{age_group_file_name}"
            df_unique.to_csv(unique_tss_file_path, sep='\t', index=False)
            print(f"Done saving the file {unique_tss_file_path}")
            print(f"\n_____________________________\n")


merging()



