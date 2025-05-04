import numpy as np
from sklearn.mixture import GaussianMixture
from .gmm_fixed_means import GaussianMixtureFixedMeans

class GaussianMixtureModel:
    def __init__(self, n_components=1, covariance_type='full'):
        
        self.model = GaussianMixture(n_components=n_components, covariance_type=covariance_type)

    def fit(self, X):
        self.model.fit(X)

    def predict(self, X):
        return self.model.predict(X)

    def score(self, X):
        return self.model.score(X)

class GaussianMixtureModelFixedMeans:
    def __init__(self, n_components=1, covariance_type='full', means_init = np.array([0.5]).reshape(-1,1)):
        
        self.model = GaussianMixtureFixedMeans(n_components=n_components, covariance_type=covariance_type, means_init = means_init)

    def fit(self, X):
        self.model.fit(X)

    def predict(self, X):
        return self.model.predict(X)

    def score(self, X):
        return self.model.score(X)