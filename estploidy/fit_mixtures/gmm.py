import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

"""
Gaussian Mixture Model (GMM) using PyTorch
This code implements a simple Gaussian Mixture Model (GMM) using PyTorch.
It includes methods for fitting the model to data and calculating the likelihood of new data points.
"""

class GaussianMixture(nn.Module):
    def __init__(self, n_components):
        super(GaussianMixture, self).__init__()
        self.n_components = n_components
        
        # Initialize parameters
        self.weights = nn.Parameter(torch.ones(n_components) / n_components)
        self.means = nn.Parameter(torch.randn(n_components))
        self.stds = nn.Parameter(torch.ones(n_components))
        
    def forward(self, x):
        # Expand dimensions for broadcasting
        x = x.unsqueeze(-1)  # Shape: (n_samples, 1)
        
        # Calculate Gaussian probability for each component
        gaussian_prob = torch.exp(-0.5 * ((x - self.means) / self.stds)**2) / (self.stds * torch.sqrt(torch.tensor(2 * np.pi)))
        
        # Weight the probabilities
        weighted_prob = self.weights * gaussian_prob
        
        # Sum probabilities for all components
        return torch.sum(weighted_prob, dim=1)
    
    def fit(self, data, n_iterations=1000, lr=0.01):
        optimizer = optim.Adam(self.parameters(), lr=lr)
        x = torch.from_numpy(data).float()
        
        for i in range(n_iterations):
            optimizer.zero_grad()
            likelihood = self.forward(x)
            loss = -torch.mean(torch.log(likelihood + 1e-10))  # Negative log likelihood
            loss.backward()
            optimizer.step()
            
            # Ensure weights sum to 1
            with torch.no_grad():
                self.weights.data = torch.softmax(self.weights.data, dim=0)
                self.stds.data = torch.abs(self.stds.data)
                
            if i % 100 == 0:
                print(f'Iteration {i}, Loss: {loss.item():.4f}')