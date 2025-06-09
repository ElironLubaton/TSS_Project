import pandas as pd
import gffpandas.gffpandas as gffpd
import numpy as np

"""
Things left to do:
add a new column named "true_start":
When the strand is +, the 'true_start' will be 'start'.
When the strand is -, the 'true_start' will be 'end'.
"""

# def loading_gene_annotation():
#     """
#     This function loads the gene annotation file, and using right now only 4 annotations:
#     - transcript_id: The id of the transcript.
#     - strand:        The direction of the strand - can be either + or -
#     - start:         The transcript start site position (depends on the strand's direction).
#     - end:           The transcript end site position   (depends on the strand's direction).
#
#     *Note - There are many more annotations to use, but I'm using only these 3 for now.
#     *Note - The origin of the annotation file is GENCODE - https://www.gencodegenes.org/human/.
#     """
#
#     raw_data_path = '/content/drive/MyDrive/Bioinformatics/gencode.v48.annotation.gff3'
#     annotation = gffpd.read_gff3(raw_data_path)
#
#     # Deleting all the numbers after the decimal number in the 'gene_id' column
#     acc_n = annotation.filter_feature_of_type(['gene']).attributes_to_columns().gene_id.str.split('.', expand=True)[0]
#
#     # Choosing the columns I want
#     gene_annotation = annotation.filter_feature_of_type(['transcript']).attributes_to_columns()[
#         ['transcript_id', 'start', 'end', 'strand']].join(acc_n)
#
#     # Cleaning the data - in 'transcript_id' column, deleting the decimal number and all the numbers after it
#     gene_annotation['transcript_id'] = gene_annotation['transcript_id'].str.replace(r'\.[A-Za-z0-9]+$', '', regex=True)
#
#     # Remove the last column of df - because it's just filled with NaNs
#     gene_annotation.drop(gene_annotation.columns[-1], axis=1, inplace=True)
#
#
#     ##  Shaked's original code:
#     # raw_data_path = '/Users/shakedshanas/Desktop/HaifaUni/Masters/div_eveness/GTEx/rawGTExdata/'
#     # annotation = gffpd.read_gff3(raw_data_path + 'gencode.v26.annotation.gff3')
#     # acc_n = annotation.filter_feature_of_type(['gene']).attributes_to_columns().gene_id.str.split('.', expand=True)[0]
#     # gene_annotation = annotation.filter_feature_of_type(['gene', 'exon']).attributes_to_columns()[['seq_id', 'type', 'transcript_type', 'gene_name', 'gene_type', 'start', 'end']].join(acc_n)
#     # gene_annotation.rename(columns={0: 'accession_number'}, inplace=True)                              # Creating the accesion_number column
#     # gene_annotation = gene_annotation[~gene_annotation.seq_id.isin(['chrX', 'chrY', 'chrM'])]          # Filtering genes that are inside chromoses X, Y and Mitochond
#     # gene_annotation = gene_annotation[gene_annotation.transcript_type.isin(['protein_coding',None])]   # Filtering out genes that are NOT protein coding









