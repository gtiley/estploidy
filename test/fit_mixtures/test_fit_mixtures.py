import numpy as np
import tempfile
from estploidy.fit_mixtures.fit_mixtures import fit_gmm_to_ab

def test_fit_gmm_to_ab():
    """
    Test that primary optimization function can return expected results from some simulated data
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        ploidy_results = []
        ab_left = np.random.normal(0.25, 0.05, (500, 4))
        ab_middle = np.random.normal(0.5, 0.05, (1000, 4))
        ab_right = np.random.normal(0.75, 0.05, (500,4))
        allele_balance_array = np.concatenate([ab_left, ab_middle, ab_right], axis=0)
        #print(allele_balance_array)
        allele_mask_array = np.random.randint(0, 2, size=(2000, 4))
        ab_dat = np.array([allele_balance_array,allele_mask_array])
        #print(ab_dat)

        for i in range(len(ab_dat[0,0,:])):
            print(i)
            ind_dat = ab_dat[0,:,i]
            #print(ind_dat)
            ind_mask = (ab_dat[1,:,i] == 1)
            #print(ind_mask)
            ind_dat_filtered = ind_dat[ind_mask].reshape(-1, 1)
            best_n = fit_gmm_to_ab(ind_name = f'{i}', dat = ind_dat_filtered, ploidy = [2,3,4,5,6], model_constraints = 1, output_dir = temp_dir)
            ploidy_results.append(best_n)
        assert ploidy_results == [4,4,4,4]