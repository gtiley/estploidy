import numpy as np
from estploidy.fit_mixtures.plot_mixtures import plot_gmm_fit_sklearn
from estploidy.fit_mixtures.gmm import GaussianMixture

def fit_gmm_to_ab(dat):
    """
    Fit Gaussian Mixture Model (GMM) to allele balance data.
    
    Parameters:
        dat (np.array): Allele balance data.
        n_components (int): Number of components in the GMM.
    """
    # Fit GMM to allele balance data
    best_n = 1
    best_bic = -np.inf
    for i in range(1, 7):
        gmm = GaussianMixture(n_components= i)
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
        plot_gmm_fit_sklearn(dat, gmm, title=f"GMM Fit to Allele Balance Data (Sample {i+1})")

    return(best_n)

def est_ploidy(ab_dat, method):
    """
    Estimate ploidy from allele balance data using the specified method.
    
    Parameters:
        ab_dat (np.array): Allele balance data returned from get_ind_freqs.
        method (str): Method for estimating ploidy ('gmm' or 'other').
    
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
            best_n = fit_gmm_to_ab(dat)
            print(f'{best_n}')
            
    else:
        raise ValueError("Unsupported method. Use 'gmm'.")