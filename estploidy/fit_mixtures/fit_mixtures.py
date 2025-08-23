import numpy as np
import sys
import math
from estploidy.fit_mixtures.plot_mixtures import plot_gmm_fit_sklearn
#from estploidy.fit_mixtures.gmm import GaussianMixture
from estploidy.fit_mixtures.gmm import GaussianMixtureFixedMeans

def get_fixed_params(n_components):
    means = None
    if (n_components > 6):
        sys.ext('ERROR: Ploidy greater than 6 is not implemented or recommended! Stopping.')
    elif(n_components < 1):
        sys.exit('ERROR: The maximum ploidy must be a positive integer! Stopping.')
    else:
        if n_components == 1:
            means = np.array([0.5]).reshape(-1,1)
        if n_components == 2:
            means = np.array([1/3,2/3]).reshape(-1,1)
        if n_components == 3:
            means = np.array([0.25,0.5,0.75]).reshape(-1,1)
        if n_components == 4:
            means = np.array([0.2,0.4,0.6,0.8]).reshape(-1,1)
        if n_components == 5:
            means = np.array([1/6,2/6,0.5,4/6,5/6]).reshape(-1,1)
    return(means)

def fit_gmm_to_ab(ind_name, dat, max_ploidy, plot_name, output_dir):
    """
    Fit Gaussian Mixture Model (GMM) to allele balance data.
    
    Parameters:
        ind_name (string): The name of the individual
        dat (np.array): Allele balance data.
        n_components (int): Number of components in the GMM.
    """
    # Fit GMM to allele balance data
    best_n = 1
    best_bic = -np.inf
    best_gmm = None
    for i in range(1, max_ploidy):
        #gmm = GaussianMixture(n_components= i)
        means = get_fixed_params(i)
        gmm = GaussianMixtureFixedMeans(n_components = i, means_init = means)
        gmm.fit(dat)
        score = gmm.score(dat)
        bic = gmm.bic(dat)
        output_file = f'{output_dir}/{ind_name}.fit.txt'
        outfile = open(output_file, 'w')
        outfile.write("Fitted GMM parameters:")
        outfile.write(f'Means:\n {gmm.means_}')
        outfile.write(f'Covariances:\n {gmm.covariances_}')
        outfile.write(f'Weights:\n {gmm.weights_}')
        # Print the best likelihood
        outfile.write(f'Best likelihood: {score}')
        outfile.write(f'BIC: {bic}')
        outfile.close()
        if bic > best_bic:
            best_bic = bic
            best_n = i
            best_gmm = gmm
    plot_gmm_fit_sklearn(dat, best_gmm, output_dir, plot_name=f'{ind_name}.fit', title=f'GMM Fit to Allele Balance Data ({ind_name})')

    return(best_n)

def est_ploidy(tax_list, ab_dat, method, max_ploidy, minimum_sites, output_dir):
    """
    Estimate ploidy from allele balance data using the specified method.
    
    Parameters:
        tax_list (list): A list of individual names corresponding to the individual order of ab_dat
        ab_dat (np.array): Allele balance data returned from get_ind_freqs.
        method (str): Method for estimating ploidy ('gmm' or 'other').
        max_ploidy (int): The maximum ploidy expected.
    
    Returns:
        ploidy_df: DataFrame containing estimated ploidy for each individual.
    """
    if method == 'gmm':
        output_file = f'{output_dir}/ploidy.txt'
        outfile = open(output_file, 'w')
        for i in range(len(ab_dat[0,0,:])):
            ind_name = tax_list[i]
            ind_dat = ab_dat[0,:,i]
            ind_mask = ab_dat[3,:,i] == 1
            ind_dat_filtered = ind_dat[ind_mask]
            ind_dat_buffer = (ind_dat_filtered > 0.05) & (ind_dat_filtered < 0.95)
            ind_dat_filtered_truncated = ind_dat_filtered[ind_dat_buffer]
            dat = ind_dat_filtered_truncated
            if len(ind_dat_filtered_truncated.shape) == 1:
                dat = ind_dat_filtered_truncated.reshape(-1, 1)
            # Fit GMM to allele balance data
            print(f"Individual {ind_name}:")
            #print(f'min sites: {minimum_sites}\n')
            #print(f'{len(dat[:])}\t{dat.shape}\n')
            if len(dat[:]) >= minimum_sites:
                best_n = fit_gmm_to_ab(ind_name, dat, max_ploidy, i, output_dir)
                #print(f'{best_n}')
                outfile.write(f'{ind_name}\t{best_n + 1}\n')
            else:
                print('Warning: Sample skipped due to low site count passing filters.\n')
        outfile.close()
    else:
        raise ValueError("Unsupported method. Use 'gmm'.")