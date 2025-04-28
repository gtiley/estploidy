import numpy as np
import pandas as pd
from estploidy.fit_mixtures.gmm import GaussianMixture

def fit_gmm_to_ab(ind_dat, n_components):
    """
    Fit Gaussian Mixture Model (GMM) to allele balance data.
    
    Parameters:
        ab_dat (np.array): Allele balance data.
        n_components (int): Number of components in the GMM.
    """
    dat = ind_dat
    # Reshape the array to 2D if necessary
    # Example: allele_balance_array = np.random.rand(100, 10)  # Replace with actual data
    if len(ind_dat.shape) == 1:
        dat = ind_dat.reshape(-1, 1)
    # Fit GMM to allele balance data   
    gmm = GaussianMixture(n_components=n_components)
    gmm.fit(dat)
    # Print fitted parameters
    print("Fitted GMM parameters:")
    print(f"Means: {gmm.means.data.numpy()}")
    print(f"Std devs: {gmm.stds.data.numpy()}")
    print(f"Weights: {gmm.weights.data.numpy()}")
    return(gmm)

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
            # Fit GMM to allele balance data
            print(f"Individual {i}:")
            gmm = fit_gmm_to_ab(ind_dat_filtered, n_components=2)
            # Save the fitted GMM parameters
            # (e.g., means, std devs, weights) for each individual
            # You can store these in a DataFrame or any other structure as needed
            # For example, you can create a DataFrame to store the results
            # ploidy_df = pd.DataFrame({'Means': gmm.means.data.numpy(),
            #                           'Std devs': gmm.stds.data.numpy(),
            #                           'Weights': gmm.weights.data.numpy()})
            # Note: You may want to modify the above line to store results in a more structured way
            # For now, just print the results
    else:
        raise ValueError("Unsupported method. Use 'gmm'.")
