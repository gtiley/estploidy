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

def fit_gmm_to_ab(dat, max_ploidy, plot_name, output_dir):
    """
    Fit Gaussian Mixture Model (GMM) to allele balance data.
    
    Parameters:
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
        print("Fitted GMM parameters:")
        print("Means: ", gmm.means_)
        print("Covariances: ", gmm.covariances_)
        print("Weights: ", gmm.weights_)
        # Print the best likelihood
        print(f'Best likelihood: {score}')
        print(f'BIC: {bic}')
        if bic > best_bic:
            best_bic = bic
            best_n = i
            best_gmm = gmm
    plot_gmm_fit_sklearn(dat, best_gmm, output_dir, plot_name, title=f'GMM Fit to Allele Balance Data (Sample {i+1})')

    return(best_n)

def est_ploidy(ab_dat, method, max_ploidy, output_dir):
    """
    Estimate ploidy from allele balance data using the specified method.
    
    Parameters:
        ab_dat (np.array): Allele balance data returned from get_ind_freqs.
        method (str): Method for estimating ploidy ('gmm' or 'other').
        max_ploidy (int): The maximum ploidy expected.
    
    Returns:
        ploidy_df: DataFrame containing estimated ploidy for each individual.
    """
    if method == 'gmm':
        for i in range(len(ab_dat[0,0,:])):
            ind_dat = ab_dat[0,:,i]
            ind_mask = ab_dat[3,:,i] == 1
            ind_dat_filtered = ind_dat[ind_mask]
            dat = ind_dat_filtered
            if len(ind_dat_filtered.shape) == 1:
                dat = ind_dat_filtered.reshape(-1, 1)
            # Fit GMM to allele balance data
            print(f"Individual {i}:")
            best_n = fit_gmm_to_ab(dat, max_ploidy, i, output_dir)
            print(f'{best_n}')
            
    else:
        raise ValueError("Unsupported method. Use 'gmm'.")