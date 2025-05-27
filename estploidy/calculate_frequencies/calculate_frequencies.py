import re
import numpy as np
import pandas as pd


def get_ind_freqs(ind_map, vcf_file, min_depth, min_count, min_qual, pate_flag, output_dir):
    '''
    Returns an np.array object of allele balance across sites for each individual from a multisample vcf.

    Parameters:
        vcf (string): a multisample vcf file uncompressed

    Returns:
        abmat (np.array): a matrix as an np.array
    '''
    
    allele_balance_data = {}
    site_depth_data = {}
    genotype_quality_data = {}
    passing_filter_data = {}
    chromosome_data = {}
    site_position_data = {}
    vcf_map = {}
    tax_list = []
    n_tax = 0
    n_variants = {}
    skip_header = 1

    with open(vcf_file,'r') as fh:
        for line in fh:
            line = line.strip()
            if '#CHROM' in line:
                temp = line.split()
                #print(line)
                for i in range(9, len(temp)):
                    this_tax = ''
                    if '/' in temp[i]:
                        if pate_flag == False:
                            tax_path = temp[i].split('/')
                            this_tax = tax_path[-1]
                        elif pate_flag == True:
                            tax_path = temp[i].split('/')
                            this_tax = tax_path[-2]
                    else:
                        this_tax = temp[i]
                    #print(this_tax)
                    if this_tax in ind_map.keys():
                        tax_list.append(this_tax)
                        n_variants[this_tax] = 0
                        n_tax = n_tax + 1
                        vcf_map[i] = this_tax
                        allele_balance_data[this_tax] = []
                        site_depth_data[this_tax] = []
                        genotype_quality_data[this_tax] = []
                        passing_filter_data[this_tax] = []
                        chromosome_data[this_tax] = []
                        site_position_data[this_tax] = []
                        #print(f'{this_tax}')
                skip_header = 0
            else:
                if skip_header == 0:
                    temp = line.split()
                    if (temp[6] == 'PASS' or ((pate_flag == True) and temp[6] == '.')):
                        for i in range(9, len(temp)): ####Continue fixing here
                            if i in vcf_map.keys():
                                ref_counts = 0
                                alt_counts = 0
                                total_count = 0
                                total_count = 0
                                allele_balance = 0
                                genotype_quality = 0
                                indicator = 0
                                genotype_string = temp[i]
                                genotype_fields = genotype_string.split(':')
                                # Anticipating that lines with fewer than 5 fields are not biallelic snps
                                if (len(genotype_fields) >= 5):
                                    if re.match(r'\d+,\d+', genotype_fields[1]):
                                        allele_counts = genotype_fields[1].split(',')
                                        #print(allele_counts)
                                        ref_counts = int(allele_counts[0])
                                        alt_counts = int(allele_counts[1])
                                        total_count = ref_counts + alt_counts
                                        if re.match(r'\d+', genotype_fields[3]):
                                            genotype_quality = int(genotype_fields[3])
                                        if (total_count > 0):
                                            allele_balance = alt_counts / total_count
                                        if ((total_count >= min_depth) and (ref_counts >= 1) and (alt_counts >= min_count) and (genotype_quality >= min_qual)):
                                            indicator = 1
                                    else:
                                        print(f'WARNING: Incorrectly formatted VCF fields!\n--> {vcf_map[i]} at variant {n_variants[vcf_map[i]]}\n-->{temp[0]}: {temp[1]}\n')
                                
                                allele_balance_data[vcf_map[i]].append(allele_balance)
                                site_depth_data[vcf_map[i]].append(total_count)
                                genotype_quality_data[vcf_map[i]].append(genotype_quality)
                                passing_filter_data[vcf_map[i]].append(indicator)
                                chromosome_data[vcf_map[i]].append(temp[0])
                                site_position_data[vcf_map[i]].append(temp[1])
                                n_variants[vcf_map[i]] = n_variants[vcf_map[i]] + 1

    if (output_dir != 'dummy'):
        for i in range(0, len(tax_list)):
            output_file = f'{output_dir}/{tax_list[i]}.txt'
            outfile = open(output_file, 'w')
            outfile.write('chr\tpos\tallele_balance\tdepth\tgenotype_quality\tpass_filters\n')
            for j in range(0, n_variants[tax_list[i]]):
                outstring = f'{chromosome_data[tax_list[i]][j]}\t{site_position_data[tax_list[i]][j]}\t{allele_balance_data[tax_list[i]][j]}\t{site_depth_data[tax_list[i]][j]}\t{genotype_quality_data[tax_list[i]][j]}\t{passing_filter_data[tax_list[i]][j]}\n'
                outfile.write(outstring)
            outfile.close()

    allele_balance_array = np.array(list(allele_balance_data.values())).transpose()
    site_depth_array = np.array(list(site_depth_data.values())).transpose()
    genotype_quality_array = np.array(list(genotype_quality_data.values())).transpose()
    passing_filter_array = np.array(list(passing_filter_data.values())).transpose()
    ab_dat = np.array([allele_balance_array, site_depth_array, genotype_quality_array, passing_filter_array])
    ab_df = pd.DataFrame(allele_balance_data)
    print(ab_df)
    print(ab_dat)
    print(ab_dat.shape)
    return(ab_dat)


def get_pop_freqs (ind_map, vcf_file, min_depth, min_count, output_file):
    pass
