import numpy as np
import sys
import pandas as pd
import logging
from estploidy.fit_mixtures.plot_mixtures import plot_gmm_fit_sklearn
from estploidy.fit_mixtures.gmm import GaussianMixture
from estploidy.fit_mixtures.gmm import GaussianMixtureFixedMeans
from estploidy.fit_mixtures.gmm import GaussianMixtureFixedMeansFixedWeights

def get_fixed_params(n_components):
    means = None
    if (n_components > 6):
        sys.ext('ERROR: Ploidy greater than 6 is not implemented or recommended! Stopping.')
    elif(n_components < 1):
        sys.exit('ERROR: The maximum ploidy must be a positive integer! Stopping.')
    else:
        if n_components == 1:
            means = np.array([0.5]).reshape(-1,1)
            weights = np.array([1.0])
        if n_components == 2:
            means = np.array([1/3,2/3]).reshape(-1,1)
            weights = np.array([0.5,0.5])
        if n_components == 3:
            means = np.array([0.25,0.5,0.75]).reshape(-1,1)
            weights = np.array([0.25,0.5,0.25])
        if n_components == 4:
            means = np.array([0.2,0.4,0.6,0.8]).reshape(-1,1)
            weights = np.array([0.25,0.25,0.25,0.25])
        if n_components == 5:
            means = np.array([1/6,2/6,0.5,4/6,5/6]).reshape(-1,1)
            weights = np.array([1/6,1/6,2/6,1/6,1/6])
    return(means, weights)

def fit_gmm_to_ab(ind_name, dat, ploidy, model_constraints, output_dir):
    """
    Fit Gaussian Mixture Model (GMM) to allele balance data.
    
    Parameters:
        ind_name (string): The name of the individual
        dat (np.array): Allele balance data.
        n_components (int): Number of components in the GMM.
    """
    # Fit GMM to allele balance data
    best_n = 1
    best_bic = np.inf
    best_gmm = None
    output_file = f'{output_dir}/{ind_name}.fit.txt'
    outfile = open(output_file, 'w')
    for i in range(0, len(ploidy)):
        n_components = ploidy[i] - 1
        gmm = None
        means, weights = get_fixed_params(n_components)
        if model_constraints == 0:
            gmm = GaussianMixture(n_components = n_components)
        if model_constraints == 1:
            gmm = GaussianMixtureFixedMeans(n_components = n_components, means_init = means)
        if model_constraints == 2:
            gmm = GaussianMixtureFixedMeansFixedWeights(n_components = n_components, means_init = means, weights_init = weights)
        gmm.fit(dat)
        score = gmm.score(dat)
        bic = gmm.bic(dat)
        outfile.write(f'Model for ploidy = {ploidy[i]}\n')
        outfile.write("Fitted GMM parameters:\n")
        outfile.write(f'Means:\n {gmm.means_}\n')
        outfile.write(f'Covariances:\n {gmm.covariances_}\n')
        outfile.write(f'Weights:\n {gmm.weights_}\n')
        # Print the best likelihood
        outfile.write(f'Best likelihood: {score}\n')
        outfile.write(f'BIC: {bic}\n')
        outfile.write('\n')
        # only consider a 3.2 point difference via Kass and Raftery 1995
        if bic < (best_bic - 3.2):
            best_bic = bic
            best_n = ploidy[i]
            best_gmm = gmm
        #We can return the categories for each point based on posterior probabilities too
        #Will be used in downstream linear models
        #Create permutation test to check if model is actually a good fit
        #predictions = gmm.predict(dat)
        #print(predictions)
    outfile.close()
    plot_gmm_fit_sklearn(dat, best_gmm, output_dir, plot_name=f'{ind_name}.fit', title=f'GMM Fit to Allele Balance Data ({ind_name})')

    return(best_n)

def est_ploidy(tax_list, ab_dat, method, ploidy_levels, minimum_sites, model_constraints, output_dir):
    """
    Estimate ploidy from allele balance data using the specified method.
    
    Parameters:
        tax_list (list): A list of individual names corresponding to the individual order of ab_dat
        ab_dat (np.array): Allele balance data returned from get_ind_freqs.
        method (str): Method for estimating ploidy ('gmm' or 'other').
        ploidy_levels (str): The ploidies to test passed as a comma-separated list.
        minimum_sites (int): The minimum number of sites to be considered for analysis
        model_constraints (int): The parameters to contrain where 0 is none, 1 is means, and 2 is means and weights
        output_dir (str): The output directory where all results will be directed
    
    Returns:
        ploidy_df: DataFrame containing estimated ploidy for each individual.
    """
    if method == 'gmm':
        output_file = f'{output_dir}/ploidy.txt'
        outfile = open(output_file, 'w')
        ploidy_dict = {}
        ploidy_level_list = ploidy_levels.split(',')
        ploidy = [int(p) for p in ploidy_level_list]
        logging.info(f'Testing for ploidy with the following values:\n{ploidy}\n')
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
            #print(f'min sites: {minimum_sites}\n')
            #print(f'{len(dat[:])}\t{dat.shape}\n')
            n_sites = len(dat[:])
            if n_sites >= minimum_sites:
                logging.info(f"Individual {ind_name}: {n_sites} sites")
                best_n = fit_gmm_to_ab(ind_name, dat, ploidy, model_constraints, output_dir)
                #print(f'{best_n}')
                outfile.write(f'{ind_name}\t{best_n}\n')
                ploidy_dict[ind_name] = best_n
            else:
                logging.warning(f'Individual {ind_name}: Sample skipped due to low site count passing filters.\n')
                ploidy_dict[ind_name] = None
        outfile.close()
        ploidy_df = pd.DataFrame.from_dict(ploidy_dict, orient = 'index')
        ploidy_df.reset_index(inplace=True)
        ploidy_df.columns = ['Individual','Ploidy']
        return(ploidy_df)
    else:
        logging.error('Terminated due to unavailable estimation method!\n')
        raise ValueError("Unsupported method. Use 'gmm'.")