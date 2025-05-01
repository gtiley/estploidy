from sklearn.mixture import GaussianMixture

class GaussianMixtureModel:
    def __init__(self, n_components=1, covariance_type='full'):
        
        self.model = GaussianMixture(n_components=n_components, covariance_type=covariance_type)

    def fit(self, X):
        self.model.fit(X)

    def predict(self, X):
        return self.model.predict(X)

    def score(self, X):
        return self.model.score(X)