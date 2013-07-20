import numpy as np


class Gaussian:
    """
    A Gaussian kernel.
    """
    
    def __init__(self, sigma=0.5):
        self.sigma = sigma
    
    def __call__(self, X, Z):
        return Gaussian.kernel(X, Z, self.sigma)
    
    @classmethod
    def kernel(cls, X, Z, sigma):
        """
        Computes the Gaussian kernel for the matrices X and Z.
        """
        X, Z = np.matrix(X, dtype="float32"), np.matrix(Z, dtype="float32")
        n, m = X.shape[0], Z.shape[0]
        XX, ZZ = np.multiply(X, X).sum(axis=1), np.multiply(Z, Z).sum(axis=1)
        d = np.tile(XX, (1, m)) + np.tile(ZZ.T, (n, 1)) - 2 * X * Z.T
        return np.array(np.exp(-d.T / (2. * sigma * sigma)))
