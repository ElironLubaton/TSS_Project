import pandas as pd

"""
In this section, we filter the metadata and make a new compact .txt file named 'Sample_Attributes_Phenotypes'.
That file contains only 5 columns:
- SAMPID - The ID of the sample.
- SMTS - The tissue the sample has been taken from.
- SMTDS - The sub-tissue the sample has been taken from
- AGE - The age group the sample belongs to.
- SEX - The sex of the donor of the sample.
"""

def extract_subjid(sampid):
    """ Function to extract correct SUBJID based on SAMPID pattern"""
    # After 'GTEX-' (5 characters), check the 6th character
    after_prefix = sampid[5]
    if after_prefix.isdigit():
        return sampid[:10]  # GTEX-12345
    else:
        return sampid[:9]   # GTEX-ABCDE


def metadata_processing():
    """
    Preprocessing the metadata. The original file holds many features, but we use only 3 for now:
    - SAMPID - The ID of the sample.
    - SMTS   - The tissue the sample has been taken from.
    - SMTDS  - The sub-tissue the sample has been taken from
    """

    # Loading the GTEx SampleAttributes file with specific columns
    metadata_file_path = data_raw_path + 'SampleAttributesDS_GTEx_Analysis_v10_Annotations.txt'
    metadata_df = pd.read_csv(metadata_file_path, sep='\t', usecols=["SAMPID", "SMTS", "SMTSD"])

    # Filtering out all the samples that have BMS id because they don't have age group at the SubjectPhenotype file
    metadata_df = metadata_df[~metadata_df['SAMPID'].str.startswith('BMS')]

    # Filtering out all the samples that have K-562 id because they belong to leukemia cell line
    metadata_df = metadata_df[~metadata_df['SAMPID'].str.startswith('K-562')]

    # Adding a new feature in the metadata_df - SUBJID - holds only the exact identifier (This feature will be dropped shortly after)
    metadata_df['SUBJID'] = metadata_df['SAMPID'].apply(extract_subjid)

    return metadata_df


def subject_phenotypes_processing():
    """
    Preprocessing the subject_phenotypes data. The original file holds more features, but we use only 3 for now:
    - SUBJID - The sample's ID prefix (9-10 characters).
    - SEX    - The sample's donor sex (male of female).
    - AGE    - The sample's donor age group.
    """

    # Loading the SubjectPhenotype file with specific columns
    pheno_file_path = data_raw_path + 'SubjectPhenotypesDS_GTEx_Analysis_v10_Annotations.txt'
    subject_phenotypes_df = pd.read_csv(pheno_file_path, sep='\t', usecols=["SUBJID", "SEX", "AGE"])

    # Removing rows where the AGE column has a missing value
    subject_phenotypes_df.dropna(subset=['AGE'], inplace=True)

    return subject_phenotypes_df


def tpm_reads_processing():
    """ Preprocessing the tpm/reads data """

    # Loading ONLY the columns names - all the samples IDs
    tpm_data_path = data_raw_path + 'Transcripts_TPM_GTEx_Analysis_v10_RSEMv1.3.3.txt.gz'
    samples_ids = pd.read_csv(tpm_data_path, delimiter='\t', nrows=0)

    # dropping 'transcript_id' and 'gene_id' - we want only the samples IDs
    samples_ids = samples_ids.iloc[:, 2:]

    return samples_ids



def main_preprocess():
    metadata_df = metadata_processing()                        # processing the metadata
    subject_phenotypes_df = subject_phenotypes_processing()    # processing the subject phenotypes data
    samples_ids = tpm_reads_processing()                       # processing the TPM/READS data

    # merging the dataframes based on SUBJID
    merged_df = metadata_df.merge(subject_phenotypes_df, on='SUBJID', how='left')
    # Removing 'SUBJID' column after the mergre
    merged_df = merged_df.drop(columns=['SUBJID'])


    # Filtering out from the Sample_Attributes_Phenotypes file all the SAMPID that is NOT found on the TPM/READS file.
    filtered_df = merged_df[merged_df['SAMPID'].isin(samples_ids)]
    # Resetting the df's indexes
    filtered_df.reset_index(drop=True, inplace=True)

    return filtered_df


# # path when using colab
# project_raw_path = '/content/drive/MyDrive/Bioinformatics/'
# path when using PC
project_raw_path = 'D:/Studies/MSc/First year/2 - Second semester/Research methods in Bioinformatics/Final_Project/'

data_raw_path = project_raw_path + 'GTEx_Data/RawData_v10/'

# Processing the data
output_file = main_preprocess()

# Saving the filtered dataframe
output_path = project_raw_path + '/CodeForProject/Output_Files/'
output_file.to_csv(output_path + 'Sample_Attributes_Phenotypes.txt', sep='\t', index=False)