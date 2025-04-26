import pandas as pd
import os

def check_dir(my_dir):
    if os.path.exists(my_dir):
        print('Existing output directory found. Writing output to: {my_dir}\n')
    else:
        os.makedirs(my_dir)
        print('Output files will be written to: {my_dir}\n')

def decompress_vcf(vcf_file):
    pass

def map_individuals(sample_sheet, vcf_file):
    '''
    Returns a dict mapping individual ids in the vcf to individual names that should be used for output. Confirms individuals are present in VCF and warns if missing.

    Parameters:
        sample_sheet (string): a sample sheet mapping individual ids to individual names and population names
        vcf (string): a multisample vcf file uncompressed

    Returns:
        ind_map (dict): a dictionary where keys are ids in the vcf and values are preferred names in downstream output
    '''
    df = pd.read_csv(sample_sheet,header=0)
    print(df.head)
    ind_map = df.set_index('individual').to_dict('index')
    print(ind_map)
    return(ind_map)