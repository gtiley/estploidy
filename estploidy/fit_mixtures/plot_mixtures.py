import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

def plot_gmm_fit_sklearn(data, gmm, output_dir, plot_name="result.fit.png", n_points=1000, title="GMM Fit to Data"):
    """
    Plot histogram of observed data with fitted GMM components from sklearn.
    
    Parameters:
        data (np.array): Original data used to fit the GMM
        gmm (GaussianMixtureModel): Fitted GMM model wrapper class
        n_points (int): Number of points for plotting the GMM curves
        title (str): Plot title
    """
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Plot histogram of observed data
    plt.hist(data, bins=50, density=True, alpha=0.5, label='Observed Data')
    
    # Generate points for plotting the GMM components
    x = np.linspace(min(data), max(data), n_points).reshape(-1, 1)
    
    # Plot individual components
    for i in range(gmm.n_components):
        # Get parameters for this component
        mu = gmm.means_[i][0]
        sigma = np.sqrt(gmm.covariances_[i][0][0])
        weight = gmm.weights_[i]
        
        # Calculate component distribution
        component = weight * norm.pdf(x, mu, sigma)
        plt.plot(x, component, '--', label=f'Component {i+1}')
    
    # Plot total mixture
    log_prob = gmm.score_samples(x)
    total = np.exp(log_prob)
    plt.plot(x, total, 'r-', label='Total Mixture', linewidth=2)
    
    plt.xlabel('Allele Balance')
    plt.ylabel('Density')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f'{output_dir}/{plot_name}.png', dpi=300)
    plt.close()