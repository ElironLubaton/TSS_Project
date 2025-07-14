import gffpandas.gffpandas as gffpd
import pandas as pd

import numpy as np
import os

import math

# # raw path for PC
# raw_path = 'D:/Studies/MSc/First year/2 - Second semester/Research methods in Bioinformatics/Final_Project/GTEx_Data/'

# # raw path for laptop
raw_path = 'C:/Users/eliro/Desktop/Studies/MSc/First year/2 - Second semester/Research methods in Bioinformatics/Final_Project/GTEx_Data/'

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

    
    # gene_annotation = annotation.filter_feature_of_type(['transcript']).attributes_to_columns()[
    #     ['seq_id', 'transcript_type', 'transcript_id', 'start', 'end', 'strand']].join(acc_n)
    # gene_annotation = gene_annotation[~gene_annotation.seq_id.isin(['chrX', 'chrY', 'chrM'])]          # Filtering genes that are inside chromoses X, Y and Mitochond
    # gene_annotation = gene_annotation[gene_annotation.transcript_type.isin(['protein_coding',None])]   # Filtering out genes that are NOT protein coding
    #     # Remove the last column of df - because it's just filled with NaNs
    # gene_annotation.drop(gene_annotation.columns[-1], axis=1, inplace=True)

    # # Adding 2 new columns which will be the start/end site according to the strand (+ or -)
    # gene_annotation['true_start'] = np.where(  # Adding true_start column
    #     gene_annotation['strand'] == '+',
    #     gene_annotation['start'],
    #     gene_annotation['end']
    # )

    # # Dropping the columns strand, end and start
    # gene_annotation.drop(columns=['seq_id', 'transcript_type', 'strand', 'end', 'start'], axis=1, inplace=True)
    
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


# merging()

def shannon_eveness(sequence, spreading_percentage):#, TPM_OR_Readcount, ceil=True, min_eve=0):
    """
    Computes array's (isoforms splicing choices in a gene) diversity.

    Args:
        sequence: A list which simulates a gene, where each cell represent a different isoform.
        spreading_percentage: An int which decides the spreading percentage we want to spread.
        TPM_OR_Readcount: A string which decides whether to read TPM files or Transcript Read Counts file.
        ceil: A boolean which decides whether to use ceiling.
        min_eve: A boolean which decides whether to calculate the minimum diversity.

    Returns:
        An number which represent the diversity of isoforms in a selected gene.
    """

    # Converting the sequence to numpy array, so the operation would be faster
    sequence = np.array(sequence)

    # Calculating how many actual splice sites in total, meaning sites with readings > 0. We want that len will be > 1
    # We don't want to use spreading on genes with 1 isoforms, this will return a value with maximum diversity - we return 0 diversity instead.
    if len(sequence[sequence>0]) <= 1:
        return 0

    # Calculating the percentage of read counts we want to spread
    p = ((sum(sequence)*(spreading_percentage/100)) / len(sequence))

    # # Decides whether to use ceiling or not.
    # # If we use TRC, and we define ceil==true, and 0<p<1 (greater than 0 in case we want no spreading), then we use ceiling
    # if TPM_OR_Readcount == 'Transcript Read counts' and ceil and (0<p<1):
    #     p = 1


    # if min_eve: # If true, then we calcualte the minimum diversity
    #     """ Minimum diversity: one isoform that contains all the reads in it.
    #     for example: [1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #     we take the minimum diversity of an arary.
    #     for example: [250, 250, 250, 250, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #     and we measure the minimum diversity as it has one isoform that contains all the reads in the array
    #     like in the first example. The first element will have the sum of reads.
    #     we add the p (spreading percentage) let's say that p=1%
    #     for example: [1000+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10]
    #     """
    #     sequence_new = [sum(sequence)+p] + [p for isoform_reads in sequence[1:]]
    # else:
    #     sequence_new = [isoform_reads+p for isoform_reads in sequence] # spreading
    sequence_new = [isoform_reads+p for isoform_reads in sequence] # spreading

    possible_ss = len(sequence_new) #computing how many POSSIBLE splice sites (including sites with 0 readings)
    total_reads = sum(sequence_new) #computing how many readings in TOTAL.


    #computing shannon's eveness index
    index_value = 0
    for splice_site in sequence_new:
        x = splice_site/(total_reads)
        if x!=0:
            index_value += (x)*(math.log(x))

    return round((-1*index_value) / (math.log(possible_ss)), 3)

    # return (round(-1*index_value, 3) / math.log(possible_ss))


def compute_div():
    age_dir_path = raw_path + 'Output Files/Age Groups/NON_UniqueTSS_or_GENE/NON_UniqueTSS/'

    for age_group_file_name in os.listdir(age_dir_path):
        # Processing only files that ends with .csv
        if age_group_file_name.endswith('.csv'):
            file_path = os.path.join(age_dir_path, age_group_file_name)
            age_group_df = pd.read_csv(file_path, delimiter='\t')

            age_group_df.drop(columns=['true_start'], axis=1, inplace=True)
            grouped_df = age_group_df.groupby('gene_id').agg(lambda x: shannon_eveness(list(x), 1)).reset_index()

            age_group_df_file_path = age_dir_path + 'NON_Unique_Evenness_spreading_1/' + f"Div_spread1_{age_group_file_name}"
            grouped_df.to_csv(age_group_df_file_path, sep='\t', index=False)
            print(f"Done saving the file {age_group_df_file_path}")


compute_div()

raw_path = 'C:/Users/eliro/Desktop/Studies/MSc/First year/2 - Second semester/Research methods in Bioinformatics/Final_Project/GTEx_Data/'
def sampid_tissue_diff():

    att_path = raw_path + f'Output Files/'
    att_df = pd.read_csv(att_path + 'Sample_Attributes_Phenotypes.txt', delimiter='\t')
    att_df.drop(columns=['SMTSD', 'SEX', 'AGE'], axis=1, inplace=True)

    tissues = ['Blood', 'Brain', 'Adipose Tissue', 'Muscle', 'Blood Vessel',
               'Heart', 'Thyroid', 'Kidney', 'Uterus', 'Vagina', 'Breast', 'Skin',
               'Salivary Gland', 'Adrenal Gland', 'Lung', 'Spleen', 'Pancreas',
               'Esophagus', 'Stomach', 'Colon', 'Small Intestine', 'Prostate',
               'Testis', 'Nerve', 'Liver', 'Pituitary', 'Ovary', 'Bladder',
               'Cervix Uteri', 'Fallopian Tube']

    for tissue in tissues:
        att_df_tissue = att_df[att_df['SMTS'] == tissue].reset_index()
        att_df_tissue.drop(columns=['SMTS'], axis=1, inplace=True)


        output_path = att_path + f"Sample_ID_tissue/{tissue}_SAMPID.csv"
        att_df_tissue.to_csv(output_path, sep='\t', index=False)


# sampid_tissue_diff()


def filt_genes_via_tss():
    """
    This function is used for filtering OUT genes that have many isoforms, but only a single TSS.
    """

    age_dir_path = raw_path + 'Output Files/Age Groups/NON_UniqueTSS_or_GENE/'

    for age_group_file_name in os.listdir(age_dir_path):
        # Processing only files that ends with .csv
        if age_group_file_name.endswith('.csv'):
            file_path = os.path.join(age_dir_path, age_group_file_name)

            age_group_df = pd.read_csv(file_path, delimiter='\t')

            # Keep only rows where 'gene_id' is duplicated (i.e., appears more than once)
            non_unique_gene_ids_df = age_group_df[age_group_df['gene_id'].duplicated(keep=False)].copy()
            non_unique_gene_ids_df.reset_index(drop=True, inplace=True)


            # Saving the NON uniques after summing
            non_unique_tss_file_path = age_dir_path + 'NON_UniqueTSS/' + f"NON_UniqueTSS_{age_group_file_name}"
            non_unique_gene_ids_df.to_csv(non_unique_tss_file_path, sep='\t', index=False)
            print(f"Done saving the file {non_unique_tss_file_path}")


filt_genes_via_tss()
