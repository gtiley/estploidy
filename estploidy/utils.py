def map_individuals(sample_sheet, vcf_file):
    '''
    Returns a dict mapping individual ids in the vcf to individual names that should be used for output. Confirms individuals are present in VCF and warns if missing.

    Parameters:
        sample_sheet (string): a sample sheet mapping individual ids to individual names and population names
        vcf (string): a multisample vcf file uncompressed

    Returns:
        ind_map (dict): a dictionary where keys are ids in the vcf and values are preferred names in downstream output
    '''

    import pandas as pd
    df = pd.read_csv(sample_sheet,header=0)
    print(df.head)

    