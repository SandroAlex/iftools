"""
===============================================================================
Basic functions written in python programming language.
===============================================================================
"""

# Load packages.
import numpy as np
import ctypes

from numpy.ctypeslib import ndpointer

# Library for C functions. Shared object.
libCfile = "/home/alex/Dropbox/repositories/iftools/iftools/base_aux.so"
libC = ctypes.CDLL(libCfile)

# Functions.
###############################################################################
def ts_clean(
        ts:np.array 
    ) -> np.array:
    """
    Delete nans in a single data array. Time series `ts` is an one-dimensional 
    numpy array.
    """

    return ts[~np.isnan(ts)]

###############################################################################
def ts_clean_pair(
        ts1:np.array,
        ts2:np.array, 
    ) -> tuple:
    """
    Delete nans in a pair of data arrays. Time series `ts1` a nd `ts2` are 
    both one-dimensional numpy arrays and they have the same size. After 
    masking out missing data, both time series can be compared to each other 
    because their new sizes are equal. 
    """

    # Mask out missing data.
    mask1 = ~np.isnan(ts1)
    mask2 = ~np.isnan(ts2)
    mask = mask1 & mask2

    return ts1[mask], ts2[mask]

###############################################################################
def bins(
        ts:np.array,
        rule:str="freedman_diaconis" 
    ) -> int:
    """
    Rule for number of bins when discretezing data. Options for `rule` 
    parameter are: (1) "basic", or (2) "freedman_diaconis".

    Reference:
    https://stadata.stackexchange.com/questions/179674/number-of-bins-when-computing-mutual-information
    """
    
    # Without not a numbers.
    tsc = ts_clean(ts)

    # Number of data points.
    n = tsc.shape[0]

    # In this case on average for two uniformly distributed random variables 
    # you will have at least 5 points for each cell of the histogram.
    if rule == "basic":
        nbins = np.int(np.floor(np.sqrt(n / 5))) 

    # There is no assumption on the distribution.
    elif rule == "freedman_diaconis":
        
        # Percentiles.
        q25 = np.percentile(tsc, 25)
        q75 = np.percentile(tsc, 75)
        IQR = q75 - q25

        # Variation between extremes.
        MAX = np.max(tsc)
        MIN = np.min(tsc)
        delta = MAX - MIN

        # Calculate rule only if IQR is not zero.
        if np.isclose(IQR, 0):
            nbins = 1

        # Rule. 
        else:
            nbins = int(np.abs(delta) / (2 * IQR * n ** (- 1.0 / 3.0)))

        # From the above `else`, nbins may be zero. In that case we put it as one.
        # It is like the time series is constant.    
        if nbins == 0:    
            nbins = 1
    
    return nbins 

###############################################################################
def probabilities(
        ts1:np.array, 
        ts2:np.array,
        rule:str="freedman_diaconis"
    ) -> tuple:
    """
    Given two time series, this function calculates estimates for marginal and
    joint probabilities.
    """

    # Without not a numbers.
    tsc1, tsc2 = ts_clean_pair(ts1, ts2)  

    # Calculate the number of bins for all histograms.
    nbins1 = bins(tsc1, rule=rule)
    nbins2 = bins(tsc2, rule=rule)
    nbins = [nbins1, nbins2]
    
    # Histograms. 1D and 2D numpy arrays.
    h1 = np.histogram(tsc1, bins=nbins1)[0]
    h2 = np.histogram(tsc2, bins=nbins2)[0]
    h12 = np.histogram2d(tsc1, tsc2, bins=nbins)[0]

    # Normalized probabilities.
    p1 = h1 / np.sum(h1)
    p2 = h2 / np.sum(h2)
    p12 = h12 / np.sum(h12)

    return p1, p2, p12    

###############################################################################
def shannon_entropy(
        p:np.array
    ) -> float:
    """
    Calculate Shannon entropy.
    """

    # Eliminate zeros. This means using the assumption 0*log2(0)=0.
    p = p[np.nonzero(p)]

    # Shannon entropy (in bit).
    h = - np.sum(p * np.log2(p))

    return h     

###############################################################################
def mutual_information(
        ts1:np.array, 
        ts2:np.array
    ) -> float:
    """
    Calculate the mutual information for the pair of time series.
    """

    # Without not a numbers.
    tsc1, tsc2 = ts_clean_pair(ts1, ts2)  

    # Marginal and joint probabilities.
    p1, p2, p12 = probabilities(tsc1, tsc2)

    # Entropy.
    h1 = shannon_entropy(p1)
    h2 = shannon_entropy(p2)
    h12 = shannon_entropy(p12)

    # Mutual information.
    mi = h1 + h2 - h12
    
    return mi

###############################################################################
def transfer_entropy(
        ts1:np.array, 
        ts2:np.array, 
        nbins:list=None,
        rule:str="freedman_diaconis",
        cboost=True
    ) -> float:        
    """
    Calculate the transfer entropy from the second time series to the first
    one (`ts2` -> `ts1`).
    """

    # Without not a numbers.
    tsc1, tsc2 = ts_clean_pair(ts1, ts2)      

    # X in the future.
    xnn = tsc1[1:]

    # X and Y in present.
    xn = tsc1[0:-1]
    yn = tsc2[0:-1]    

    # Calculate the number of bins for all histograms if these bins were not 
    # previously prescribed. Use predefined rule. 
    if nbins is None:
        nbins_xnn = bins(xnn, rule=rule)
        nbins_xn = bins(xn, rule=rule)
        nbins_yn = bins(yn, rule=rule)
        nbins = [nbins_xnn, nbins_xn, nbins_yn]    

    # Calculate the joint probability p(x_(n+1), x_n, y_n).
    data = np.column_stack((xnn, xn, yn))
    p_xnn_xn_yn = np.histogramdd(data, bins=nbins)[0]
    p_xnn_xn_yn = p_xnn_xn_yn / np.sum(p_xnn_xn_yn)

    # Calculate the joint probability p(x_n, y_n).
    p_xn_yn = np.histogram2d(xn, yn, bins=nbins[1:3])[0]
    p_xn_yn = p_xn_yn / np.sum(p_xn_yn)

    # Calculate the joint probability p(x_(n+1), x_n).
    p_xnn_xn = np.histogram2d(xnn, xn, bins=nbins[0:2])[0]
    p_xnn_xn = p_xnn_xn / np.sum(p_xnn_xn)

    # Call C ancillary function.
    if cboost:

        p_xnn_xn_yn = p_xnn_xn_yn.astype(np.float32)
        p_xnn_xn = p_xnn_xn.astype(np.float32)
        p_xn_yn = p_xn_yn.astype(np.float32)
        p_xn = p_xn.astype(np.float32)

        func = libC.transfer_entropy_summation
        func.restype = ctypes.c_float
        func.argtypes = [
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ndpointer(ctypes.c_float),
            ndpointer(ctypes.c_float),
            ndpointer(ctypes.c_float),
            ndpointer(ctypes.c_float)
        ]

        # Transfer entropy.
        te = func(
            nbins_xnn, nbins_xn, nbins_yn, 
            p_xnn_xn_yn, p_xnn_xn, p_xn_yn, p_xn
        )

    # Call Python function:
    else:

        # Summation.
        te = 0
        for i in range(nbins[0]):         # xnn.     
            for j in range(nbins[1]):     # xn.
                for k in range(nbins[2]): # yn.

                    # We are going to add only non zero terms.
                    # Assumptions:  0*log2(a/0)=0.
                    #               0*log2(0) = 0.
                    # Here we need to further investigate these 
                    # assumptions.
                    if(
                        not np.isclose(p_xnn_xn_yn[i, j, k], 0) and 
                        not np.isclose(p_xn[j], 0) and
                        not np.isclose(p_xn_yn[j, k], 0) and
                        not np.isclose(p_xnn_xn[i, j], 0)    
                    ):

                        # Argument for log2 function.
                        arg = (
                            p_xnn_xn_yn[i, j, k] * p_xn[j] /
                            p_xn_yn[j, k] / p_xnn_xn[i, j]
                        ) 
                    
                        # Accumulate for transfer entropy.
                        te = te + p_xnn_xn_yn[i, j, k] * np.log2(arg)

    return te        