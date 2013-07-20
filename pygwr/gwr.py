"""
pygwr
"""

import numpy as np
import os, sys, json
import gwstatsmodels.api as sm
import scipy.sparse
from scipy import linalg

from gaussian import Gaussian


def read_csv(path, header=True, sep="\t"):
    f = open(path)
    data = []
    header_read = False
    for l in f:
        if header and not header_read:
            header_read = True
            h = l.strip().split(sep)
            continue
        data.append(l.strip().split(sep))
    f.close()
    return h,data



class GWR:
    """
    An implementation of Geographically Weighted Regression.
    """
    def __init__(self, targets, samples, locations, family=None, kernel=None):
        """
        Initializes the GWR class using the target values, the samples, the
        locations of the target values and samples.
        """
        self.targets = targets
        # We add a column with 1's at the end for the intercept
        self.samples = np.column_stack((samples, np.ones(len(samples))))
        self.locations = locations
        self.kernel = kernel
        if family == 'poisson':
            self.family = sm.families.Poisson()
        else:
            self.family = None
    
    def estimate(self, estimation_pt, sample, bandwidth):
        """
        Performs a GWR estimation at a given point.
        estimation_pt is a list, tuple or numpy array containing the 
            coordinates of the regression point.
        sample is a list or array with a sample to compute a prediction at 
            the same time, by applying the regression parameters found by the
            regression to this value.
        bandwidth is a the kernel bandwidth.
        Returns a GWRResult object.
        """
        # Compute the weights of all known samples using the Euclidean
        # distance between the points, and a Gaussian kernel.
        if self.kernel:
            w = self.kernel(self.locations, estimation_pt, bandwidth)[0]
        else:
            w = Gaussian.kernel(self.locations, estimation_pt, bandwidth)[0]
        # Perform the weighted regression using Maximum likelihood
        if self.family is None: # Gaussian GWR
            res = sm.WLS(self.targets, self.samples, w)
            fres = res.fit()
        else: # Poisson GWR
            res = sm.GLM(endog=self.targets, exog=self.samples, family=self.family)
            fres = res.fit(data_weights=w)
        # Pack everything into a GWRResult object
        # A GWRResult allows convenient inspection of the results.
        gwr_result = GWRResult(estimation_pt, sample, fres, self)
        del(fres)
        return gwr_result
    
    def estimate_at_target_locations(self, bandwidth):
        """
        Performs a GWR estimation at all the target locations.
        """
        res = GWRResultSet()
        for i in range(len(self.locations)):
            estimation_pt = self.locations[i]
            sample = self.samples[i,:-1]
            try:
                gwr_result = self.estimate(estimation_pt, sample, bandwidth)
                gwr_result.idx = i
            except:
                gwr_result = None
            res.append(gwr_result)
        return res
    
    def global_regression(self):
        """
        Performs a global regression using the targets and samples only,
        by ignoring the locations. OLS is used for this regression.
        """
        if self.family is None:
            return sm.OLS(self.targets, self.samples).fit()
        else:
            return sm.GLM(endog=self.targets, exog=self.samples, family=self.family).fit()
        
    def aicc(self, bandwidth):
        """
        Computes the AICc criterion for the given bandwidth.
        """
        res = []
        n = len(self.locations)
        s = np.zeros(n, dtype=np.float64)
        x = np.matrix(self.samples, dtype=np.float64)
        dev = 0.0       # deviance
        # Make estimations and compute the diagonal of the hat matrix S
        for i in range(n):
            estimation_pt = self.locations[i]
            sample = self.samples[i,:-1]
            r = self.estimate(estimation_pt, sample, bandwidth)
            wi = scipy.sparse.spdiags(r.weights, 0, n, n)
            s[i] = (x[i] * np.matrix(linalg.pinv(np.array(x.T * wi * x))) * x.T * wi)[0,i]
            res.append(r.prediction)
            yhat = np.exp(r.prediction)
            y = self.targets[i]
            if y != 0: dev += 2 * y * (np.log(y) - np.log(yhat))
            dev -= 2 * (y - yhat)
            del r
        # Compute the residual sum of squares
        rss = 0.0
        sigma2_hat = 0.0
        sigma2_hat_b = 0.0
        for i in range(n):
            prediction = res[i]
            target_value = self.targets[i]
            rss += (prediction - target_value) * (prediction - target_value)
        # We don't need the GWR result anymore
        del(res)
        # Compute the sigma (st deviation of the error, see page 61 GWR book)
        sigma = np.sqrt(rss / n)
        # Finally the AICc
        aic, aicc, bic, sbic, K = 0., 0., 0., 0., 0.
        if self.family is not None:     # Poisson model deviance
            K = np.sum(s)       # effective number of parameters
            aic = dev + 2*K
            aicc = aic + 2 * ((K*(K+1)) / (n - K - 1))
            bic = dev + K*np.log(n)
        else:
            aicc = 2 * n * np.log(sigma) + n * np.log(2 * np.pi) + n * ( (n + np.sum(s)) / (n - 2 - np.sum(s)) )
        return {'aicc': aicc, 'aic': aic, 'bic': bic, 'dev': dev, 'K': K}


class GWRResult:
    """
    Represents the result of a GWR estimation at one estimation point.
    """
    def __init__(self, estimation_pt, sample, regression_result, gwr):
        self.gwr = gwr
        self.idx = None      # Index of the estimation point in GWR
        self.estimation_pt = estimation_pt
        self.sample = sample
        self.params = regression_result.params
        self.weights = regression_result.model.weights
        self.prediction = self.predict()
        self.result = regression_result
    
    def predict(self):
        """
        Predicts the regression value at the estimation point.
        """
        # If there is no sample, just return None
        if self.sample == None:
            return None
        # Add the intercept to the sample
        s = np.hstack((self.sample, np.ones(1)))
        # Compute the prediction
        return np.sum(s * self.params)
    
    



class GWRResultSet(list):
    """
    A collection of GWR results.
    """
    
    def export_csv(self, path):
        """
        Exports the GWR result set into a CSV text file.
        """
        f = open(path, "w")
        p0 = self.__getitem__(0)
        # Add the estimation point coordinates
        ndims = len(p0.estimation_pt)
        # Write out the header line.
        for i in range(ndims):
            f.write("EST_PT_%d\t" % (i+1))
        n = len(p0.params)
        # Add the sample
        for i in range(n-1):
            f.write("SAMPLE_%d\t" % (i+1))
        # Add the parameters
        for i in range(n):
            f.write("PARAM_%d\t" % (i+1))
        # Add the t-values
        for i in range(n):
            f.write("TVAL_%d\t" % (i+1))
        # Add the p-values
        for i in range(n):
            f.write("PVAL_%d\t" % (i+1))
        # Add the target and prediction
        f.write("TARGET\t")
        f.write("PREDICT\n")
        # Go through all the estimation points and create a feature for each
        for r in self:
            for i in range(ndims):
                f.write("%f\t" % r.estimation_pt[i])
            for i in range(n-1):
                f.write("%f\t" % r.sample[i])
            for i in range(n):
                f.write("%f\t" % r.params[i])
            for i in range(n):
                f.write("%f\t" % r.result.tvalues[i])
            for i in range(n):
                f.write("%f\t" % r.result.pvalues[i])
            if r.idx != None:
                f.write("%f\t" % r.gwr.targets[r.idx])
            else:
                f.write("\t")
            f.write("%f\n" % r.prediction)
        f.close()
    





