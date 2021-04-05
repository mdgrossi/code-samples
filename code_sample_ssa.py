#!/usr/bin/env python

# =========================================================================== #
# Author:    mgrossi
# Created:   8 July 2020
# 
#
# NOTES FOR CSS AND NOAA TEAM:
#
# This script was designed to perform multivariate singular spectrum analysis
# (MSSA) on hundreds of drifter velocity time series. The original idea was to
# decompose the zonal and meridional velocity time series for each drifter,
# identify the statistical signals to which I fit deterministic models (in the
# form of cosine regressions), and separate out the noise for a neural network
# to learn.
#
# As noted in the SSA docstring below, I started with a simple example of a
# univariate SSA presented on kaggle and expanded it to suit my needs. To start
# with, I had to create a multivariate version because velocity components 
# are not independent and could not handled separately. I also introduced
# several methods for pattern identification and reconstruction. If one chooses
# too liberal parameters for SSA, one can end up over-decomposing the time 
# series leading to multiple modes that really should be one single mode. In
# most applications, only one (or at most a handful) of time series are being
# analyzed at a time, so one could play around with these parameters manually.
# However, I could not do that for hundreds of drifters, so instead, I opted to
# prefer over-decomposing and then implement ways of objectively regrouping
# modes afterwards. The find_patterns() and the group() methods do just that
# (where a pattern is defined as a group of statically correlated modes.)
#
# I also added some plotting methods for visualization purposes. This routine
# originally called in a shell script, but I also have an interactive jupyter 
# notebook version in which I did a lot of experimenting, plotting, etc.
#
# =========================================================================== #

# To execute in Terminal:
# /anaconda3/bin/python mssa_extract_patterns.py

# =========================================================================== #
# Modules
import numpy as np
import pandas as pd
import os
import argparse

# =========================================================================== #
# Functions and classes
def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Function control parameters.',
                        epilog='Output in ~/.../DrifterPrediction/results/ssa/',
                        prog='mssa_extract_patterns',
                        usage='%(prog)s [arguments]')
    parser.add_argument('file', metavar='file', type=str, nargs='+',
                        help='csv file to operate on')
    parser.add_argument('-d', '--days', metavar='days', type=int, default=None,
                        help='Number of days in the time series to decompose')
    parser.add_argument('-l', '--lag', metavar='lag_window', type=float,
                        help='Lag window for composing trajectory matrix.',
                        default=0.35)
    parser.add_argument('-s', '--signal', metavar='threshold', type=float,
                      help='Fraction of variance (0,1] to be considered signal',
                      default=0.9)
    parser.add_argument('-c', '--correlation', metavar='correlation',
                      help='Minimum correlation coefficient for mode groupings',
                      type=float, default=0.5)
    parser.add_argument('-n', '--normalize', action='store_true',
                      help='Normalize data before analysis')
    parser.add_argument('-w', '--writeout', action='store_true',
                      help='Write data to file instead of printing to terminal')
    return parser.parse_args()

class SSA(object):
    
    __supported_types = (pd.Series, np.ndarray, list)
    
    def __init__(self, tseries, window, standardize=False, save_mem=True, name=None):
        """
        Decomposes the given univariate time series with a singular-spectrum
        analysis. Assumes the values of the time series are recorded at equal
        intervals.
        
        Parameters
        ----------
        tseries : Original time series, in the form of a Pandas Series, NumPy
            array, or list. 
        window : Window length. Must be an integer 2 <= window <= N/2, where N
            is the length of the time series.
        standardize : Standardize time series prior to decomposing by 
            subtracting the mean and dividing be the standard deviation.
            For multivariate time series, this is performed on each time series
            separately. If True, reconstructed time series will be returned
            un-standardized. Defaults to False.
        save_mem : Conserve memory by not retaining the elementary matrices.
            Recommended for long time series with thousands of values. Defaults
            to True.
        name : String or int to call instance. Defaults to None but will look 
            for a name attribute in tseries.
        
        Note: Even if a NumPy array or list is used for the initial time
        series, all time series returned will be in the form of a Pandas Series
        or DataFrame object.
        
        Adopted from kaggle, retrieved 25 May 2020:
        www.kaggle.com/jdarcy/introducing-ssa-for-time-series-decomposition
        """
        
        # Type-check the initial time series
        assert len(tseries.shape) == 1, \
                'Data appears to contain multiple time series. '\
                'Use MSSA() for multivariate SSA or check data format.'
        if not isinstance(tseries, self.__supported_types):
            raise TypeError('Unsupported time series object. '\
                            'Try Pandas Series, NumPy array, or list.')
        
        # Checks to save us from ourselves
        self.N = len(tseries)
        if not 2 <= window <= self.N/2:
            raise ValueError('The window length must be in the interval '\
                             '[2, N/2].')
        
        # Store data and parameters
        self.L = window
        self.orig_TS = pd.Series(tseries)
        self.K = self.N - self.L + 1
        self.name = name
        if not name:
            try:
                self.name = tseries.name
            except AttributeError:
                pass
        self.standardize = standardize
        
        # Standardize each time series
        if self.standardize:
            self.std_TS = self.standardized()
        else:
            self.std_TS = self.orig_TS
             
        # Embed the time series in a trajectory matrix
        self.X = np.array([self.std_TS.values[i:self.L+i] 
                           for i in range(0, self.K)]).T

        # Decompose trajectory matrix using Singular Value Decomposition (SVD)
        self.U, self.Sigma, VT = np.linalg.svd(self.X)
        self.d = np.linalg.matrix_rank(self.X)

        # Reconstructed time series matrix (to fill)
        self.TS_comps = np.zeros((self.N, self.d))

        if not save_mem:
            # Construct and save all the elementary matrices
            self.X_elem = np.array([self.Sigma[i]*np.outer(self.U[:,i], VT[i,:])
                                    for i in range(self.d)])

            # Diagonally average the elementary matrices, store them as columns
            # in array.
            for i in range(self.d):
                X_rev = self.X_elem[i, ::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in
                                      range(-X_rev.shape[0]+1, X_rev.shape[1])]

            # Save the (transposed) V array
            self.V = VT.T

        else:
            # Reconstruct the elementary matrices without storing them
            for i in range(self.d):
                X_elem = self.Sigma[i] * np.outer(self.U[:,i], VT[i,:])
                X_rev = X_elem[::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in
                                      range(-X_rev.shape[0]+1, X_rev.shape[1])]

            self.X_elem = 'Re-run with save_mem=False to retain the '\
                          'elementary  matrices.'

            # The V array may also be very large under these circumstances,
            # so discard it.
            self.V = 'Re-run with save_mem=False to retain the V matrix.'
                
        # Calculate the w-correlation matrix.
        self.calc_wcorr()
    
    def calc_msre(self, indices):
        """
        Calculate mean squared reconstruction error.
        """
        if isinstance(indices, int): indices = [indices]

        self.msre = (np.abs(self.orig_TS - self.reconstruct(indices))**2).mean()
        return self.msre

    def standardized(self):
        """Standardize the original time series by subtracting mean and 
        dividing by standard deviation."""
        return (self.orig_TS - self.orig_TS.mean()) / self.orig_TS.std()

    def unstandardize(self, tseries):
        """Unstandardize the output time series."""
        return (tseries * np.array(np.std(self.orig_TS, ddof=1))) + \
                np.array(np.mean(self.orig_TS))
    
    def reconstruct(self, indices):
        """
        Reconstructs the time series from its elementary components, using the
        given indices. Returns a Pandas Series object with the reconstructed
        time series.
        
        Parameters
        ----------
        indices : An integer, list of integers or slice(n,m) object,
            representing the elementary components to sum.
        """
        # Check types and convert if needed
        if isinstance(indices, int): indices = [indices]
        if isinstance(indices, tuple): indices = list(indices)
        
        # Extract and sum up the desired components
        ts_vals = self.TS_comps[:,indices].sum(axis=1)
        ts_vals = pd.DataFrame(ts_vals, index=self.orig_TS.index)
        
        # Unstandardize if necessary
        if self.standardize:
            ts_vals = self.unstandardize(ts_vals)
        
        return ts_vals
    
    def reconstruct_pattern(self, pattern):
        """
        Reconstruct the time series from the elementary components of specific 
        pattern from find_patterns().
        """
        return self.reconstruct(self.get_patterns()[pattern])
    
    def calc_wcorr(self, w_repeats=1):
        """
        Calculates the w-correlation matrix for the time series and cumulative
        weighted correlations, C(k), between all possible signal-noise
        reconstructions. See docstrings for plot_wcorr_list() and
        get_wcorr_list() for further explanations.
        """
        
        # Calculate the weights
        L_star = min(self.L, self.K)
        K_star = max(self.L, self.K)
        w = np.array(list(np.arange(L_star)+1) + \
                     [L_star]*(K_star-L_star-1) + \
                     list(np.arange(L_star)+1)[::-1])
        w = np.tile(w, w_repeats)
        
        def w_inner(F_i, F_j):
            return w.dot(F_i*F_j)
        
        # Calculate weighted norms, ||F_i||_w, then invert.
        F_wnorms = np.array([w_inner(self.TS_comps[:,i], self.TS_comps[:,i])
                             for i in range(self.d)])
        F_wnorms = F_wnorms**-0.5
        
        # Calculate Wcorr matrix.
        self.Wcorr = np.identity(self.d)
        for i in range(self.d):
            for j in range(i+1,self.d):
                self.Wcorr[i,j] = abs(w_inner(self.TS_comps[:,i],
                                self.TS_comps[:,j]) * F_wnorms[i] * F_wnorms[j])
                self.Wcorr[j,i] = self.Wcorr[i,j]
        
        # Calculate Wcorr for pairs of reconstructed components.
        F_wnorms_signal = [abs(w.dot(self.TS_comps[:,:i].sum(axis=1)))**-0.5 \
                           for i in range(1,self.d)]
        F_wnorms_noise = [abs(w.dot(self.TS_comps[:,i:].sum(axis=1)))**-0.5 \
                          for i in range(1,self.d)]
        corrs = [abs(w.dot(self.TS_comps[:,:i].sum(axis=1) * \
                           self.TS_comps[:,:i].sum(axis=1))) \
                 for i in range(1,self.d)]
        self.wcorr_list = np.array(corrs) * \
                          (np.array(F_wnorms_signal) * np.array(F_wnorms_noise))
    
    def plot_wcorr(self, min=None, max=None):
        """
        Plots the w-correlation matrix for the decomposed time series. Use 'min'
        and 'max' to zoom in.
        """
        if min is None:
            min = 0
        if max is None:
            max = self.d
        
        if self.Wcorr is None:
            self.calc_wcorr()
        
        fig, ax = plt.subplots(figsize=(10,10))
        im = ax.imshow(self.Wcorr, cmap='Blues')
        plt.xlabel(r'$\tilde{F}_i$')
        plt.ylabel(r'$\tilde{F}_j$')
        cb = plt.colorbar(im, fraction=0.045)
        cb.set_label('$W_{i,j}$')
        cb.mappable.set_clim(0,1)
        
        # For plotting purposes:
        if max == self.d:
            max_rnge = self.d-1
        else:
            max_rnge = max
        
        plt.xlim(min-0.5, max_rnge+0.5)
        plt.ylim(max_rnge+0.5, min-0.5)
    
    def plot_wcorr_list(self, xmin=None, xmax=None):
        """
        Plots a list of cumulative weighted correlations, C(k), between all
        possible signal-noise reconstructions. See get_wcorr_list() docstring 
        for explanation.
        
        From Patterson et al. (2011):
        'The existence of structure in the series is indicated by local minima
        and maxima in the graph of C(k) against k. For example, a typical
        pattern is a decline in the cumulative w-correlations corresponding to
        a separation of the components as k increases; the first local minimum
        indicates the first separation and subsequent local maxima suggest
        possible secondary structure.'
        
        Patterson, K., H. Hassani, S. Heravi, A. Zhigljavsky (2011)
            'Multivariate Singular Spectrum Analysis for Forecasting Revisions
            to Real-time Data, J. Appl. Stats.,
            DOI: 10.1080/02664763.2010.545371.
        """
        fig, ax = plt.subplots(figsize=(10,10))
        ax.plot(self.get_wcorr_list())
        ax.set_xlabel('k')
        ax.set_ylabel('Weighted correlation, C(k)')
        ax.set_xlim(xmin, xmax)
        plt.show()
    
    def plot_all_patterns(self):
        # Fiddle with color cycle - need more colors!
        plt.rcParams['font.size'] = 14
        plt.rcParams['image.cmap'] = 'plasma'
        plt.rcParams['axes.linewidth'] = 2
        
        fig, ax = plt.subplots(1, 1, figsize=(10,8))
        color_cycle = cycler(color=plt.get_cmap('tab20').colors)
        ax.set_prop_cycle(color_cycle)
        ax.plot(self.get_original(), alpha=1, lw=1, color='black',
                      linestyle='dashed', zorder=0)
        for key in self.patterns.keys():
            ax.plot(self.reconstruct(self.patterns[key]), lw=2, zorder=1)
        ax.set_xlabel("$t$")
        ax.set_ylabel(r"$\tilde{F}_i(t)$")
        legend = [r"$\tilde{F}_{%s}$" %(i-1)
                  for i in range(len(self.patterns.keys())+1)] + ["$F$"]
        ax.set_title("{} Patterns of Drifter {}"\
                       .format(len(self.patterns), self.name))
        ax.legend(legend, loc=(1.05,0.1))
        plt.show()
    
    def plot_pattern(self, pattern):
        # Fiddle with color cycle - need more colors!
        plt.rcParams['font.size'] = 14
        plt.rcParams['image.cmap'] = 'plasma'
        plt.rcParams['axes.linewidth'] = 2
        
        fig, ax = plt.subplots(1, 1, figsize=(10,8))
        ax.plot(self.get_original(), alpha=1, color='black', lw=1,
                linestyle='dashed', zorder=0, label='Original')
        ax.plot(self.reconstruct(self.get_patterns()[pattern]),
                label='Pattern', linestyle='dotted')
        for v in self.get_patterns()[pattern]:
            ax.plot(self.reconstruct(int(v)), label='Mode {}'.format(v),
                    zorder=1)
        ax.set_xlabel("$t$")
        ax.set_ylabel(r"$\tilde{F}_i(t)$")
        plt.legend(loc=(1.05,0.1))
        plt.show()
    
    def group(self, signal_cutoff):
        """
        Creates a dictionary of correlated SSA components, where each entry is
        a tuple of correlated components. Pairs are first sorted by correlation
        and then matched with the appropriate group based on highest
        correlation value.

        The components listed in each dictionary entry can be added together to 
        create a single statistical pattern of variability when reconstructing
        the time series. 
    
        Parameters
        ----------
        pairs : A numpy array or list (to be converted to numpy array)
            containing pairs of correlated indices from SSA. These pairs can be
            tuples, lists, or columns in an array, e.g.: (0, 1).
        corr : A numpy array or list (to be converted to numpy array)
            containing the correlation between each pair of indices. Must be
            the same length as the number of pairs passed.
        """
        # Make sure each pair has a correlation value for sorting
        assert len(self.pairings)==len(self.weights), \
               'Every pair must have a weight: len(pairs)==len(corr).'
        # Change pairs and corr to numpy array if needed
        if isinstance(self.pairings, list):
            self.pairings = np.array(self.pairings)
        if isinstance(self.weights, list):
            self.weights = np.array(self.weights)
        # Dictionary in which to store correlated SSA components
        fDict = {}
        # If no statistical groups exist, store each mode individually and 
        # return the resulting dictionary
        if any(elem is None for elem in self.weights):
            for i in self.pairings:
                fDict['F{}'.format(str(i).zfill(2))] = tuple([i])
            return fDict
        # Sort pairings by correlation so that for modes that are paired with 
        # more than one other mode, the strongest-correlated pair is taken
        # first.
        self.pairings = self.pairings[self.weights.argsort()[::-1]]
        # Check for signal modes that are not paired with any other mode, make
        # new dictionary entries for each. (This is often true of the most
        # dominant mode.) Accomplish this by subtracting two sets: one
        # containing all modes defined as signal, and the second set containing
        # all modes included in self.pairings. Any mode that exists in the first
        # set but not the second meets this criteria.
        sn_ind = sum(self.get_variance() <= signal_cutoff)
        to_add = list(set(np.arange(0, sn_ind)) -
                      set(np.unique(self.pairings)))
        # Counter for dictionary key
        counter = 0
        # if len(to_add) == 1:
        #     fDict['F00'] = fDict['F00'] + tuple(to_add)
        # First store non-paired signal modes, if any
        if len(to_add) > 0:
            for d in to_add:
                fDict['F{}'.format(str(counter).zfill(2))] = tuple([d])
                counter += 1
        # Now add the pairs by looking at both elements in the pair
        for i,j in self.pairings:
            # If both values are already in a dictionary set, do nothing
            if all(y in [x for value in fDict.values() for x in value]
                   for y in [i,j]):
                pass
            # If neither value is already in a dictionary set, make a new
            # key/set
            elif all(y not in [x for value in fDict.values() for x in value]
                     for y in [i,j]):
                fDict['F{}'.format(str(counter).zfill(2))] = tuple([i, j])
                counter += 1
            # If the first value is already in the dictionary, add the second
            # value to the set that contains the first value, and vice versa
            elif i in [x for value in fDict.values() for x in value]:
                for k, v in fDict.items():
                    if i in v:
                        fDict[k] = v + tuple([j])
            elif j in [x for value in fDict.values() for x in value]:
                for k, v in fDict.items():
                    if j in v:
                        fDict[k] = v + tuple([i])
    
        # Sort dictionary items according to original eigenvalues
        sorted_vals = sorted(list(fDict.values()), key=lambda x: x[0])
        for key, new_val in zip(fDict.keys(), sorted_vals):
            fDict[key] = new_val
        return fDict
    
    def find_patterns(self, signal_cutoff=0.9, corr_strength=0.5):
        """
        Group SSA elementary components by correlation to identify prominant 
        patterns of variability in time series.
        
        Parameters:
        signal_cutoff : Float (0,1] specifying the fraction of the variance to 
            to take as signal; the remainder is discarded as noise.
        corr_strength : Float (0,1) specifying the threshold for the magnitude 
            of the correlation coefficient, r, when grouping modes together.
            Correlation relationships can typically be categorized as follows:
                r < 0.3           very weak or no relationship
                0.3 <= r < 0.5    weak
                0.5 <= r < 0.7    moderate
                r >= 0.7          strong
        """
        # Signal-noise splitting index
        sn_ind = sum(self.get_variance() <= signal_cutoff)
        # Upper triangular correlation matrix
        tri = self.get_wcorr()[0:sn_ind, 0:sn_ind] * \
              np.tri(*self.get_wcorr()[0:sn_ind,0:sn_ind].shape).T
        # Group SSA elementary components by correlation
        self.pairings = np.argwhere((tri >= corr_strength) & (tri != 1))
        self.weights = tri[(tri >= corr_strength) & (tri != 1)]
        # If no statistical pairings exist, store a list of mode numbers
        if len(self.pairings) == 0:
            self.pairings = [i for i in range(sn_ind)]
            self.weights = [None] * sn_ind
        # Generate patterns
        if len(self.pairings) == 0:
            self.pairings = 'No patterns exist within specified signal cutoff.'
            print('Drifter ' + self.name + ': ', self.pairings)
            self.patterns = []
        else:
            self.patterns = self.group(signal_cutoff=signal_cutoff)
        # Return results
        return [self.get_window(), len(self.get_original()), self.K, sn_ind, 
                len(self.patterns)]
    
    def get_wcorr(self):
        """
        Returns the weighted-correlation matrix.
        """
        return self.Wcorr
            
    def get_wcorr_list(self):
        """
        Returns a list of cumulative weighted correlations between all possible
        signal-noise reconstructions, where each of the k values (1 <= k <= d-1)
        compares the two time series resulting from considering only the first
        k elementary component as signal and the rest to be noise.
        
        Adapted from:
        Patterson, K., H. Hassani, S. Heravi, A. Zhigljavsky (2011)
            'Multivariate Singular Spectrum Analysis for Forecasting Revisions
            to Real-time Data, J. Appl. Stats.,
            DOI: 10.1080/02664763.2010.545371.
        """
        return self.wcorr_list

    def get_original(self):
        """
        Returns the original time series as a Pandas series object.
        """
        return self.orig_TS
        
    def get_components(self, n=0):
        """
        Returns all the time series components in a single Pandas DataFrame
        object.
        """
        if n > 0:
            n = min(n, self.d)
        else:
            n = self.d
        
        # Create list of columns - call them F0, F1, F2, ...
        cols = ['F{}'.format(i) for i in range(n)]
        return pd.DataFrame(self.TS_comps[:, :n], columns=cols,
                            index=self.orig_TS.index)
    
    def get_window(self):
        """
        Returns the lag window used for decomposition.
        """
        return self.L

    def get_msre(self):
        """Returns the mean squared reconstruction error.
        """
        try:
            return self.msre
        except AttributeError:
            print('Nothing to return. Run calc_msre() to calculate first.')
    
    def get_variance(self, head=None):
        """Returns the cumulative total variation explained by the elementary
        matrices.
        """
        return ((self.Sigma**2).cumsum() / (self.Sigma**2).sum())[0:head]
    
    def get_patterns(self):
        """Returns dictionary of patterns"""
        try:
            return(self.patterns)
        except AttributeError:
            print('Nothing to return. Run find_patterns() first.')
        
    def __str__(self):
        return 'SSA object containing singular-spectrum-analysis results: {}'\
               .format(', '.join(list(self.__dict__.keys())))
    
    def __repr__(self):
        return '{lag_window:' + str(self.L) + ', tseries_length:' + \
                str(self.N) + ', trajMatLen:' + str(self.K)+'}'
 
class MSSA(SSA):
    
    __supported_types = (pd.DataFrame, np.ndarray)
    
    def __init__(self, tseries, window, save_mem=True, name=None):
        """
        Decomposes the given set of time series set with a multivariate
        singular-spectrum analysis. Assumes the values of the time series are
        recorded at equal intervals and are appended to each other column-wise.
        
        Parameters
        ----------
        tseries : Original set time series, in the form of a Pandas DataFrame 
            or NumPy array, with samples arranged in rows and time series in 
            columns. 
        window : Window length. Must be an integer 2 <= window <= N/2, where N
            is the length of the time series.
        save_mem : Conserve memory by not retaining the elementary matrices.
            Recommended for long time series with thousands of values. Defaults
            to True.
        name : String or int to call instance. Defaults to None but will look 
            for a name attribute in tseries.
        
        Note: Even if a NumPy array is used for the initial time series, all
        results will be returned as a Pandas Series or DataFrame object.
        
        Adopted from kaggle, retrieved 25 May 2020:
        www.kaggle.com/jdarcy/introducing-ssa-for-time-series-decomposition
        """
        
        # Type-check the initial time series
        assert len(tseries.shape) > 1, \
                'Data appears to contain a single time series. '\
                'Use SSA() for univariate SSA or check data format.'
        if not isinstance(tseries, self.__supported_types):
            raise TypeError('Unsupported time series object. '\
                            'Try Pandas Series, NumPy array or list.')
             
        # Checks to save us from ourselves
        if not 2 <= window <= tseries.size//2:
            raise ValueError('The window length must be in the interval '\
                             '[2, N/2].')
        
        # Store data and parameters
        self.N = len(tseries)
        self.L = window
        self.orig_TS = pd.DataFrame(tseries)
        self.K = self.N - self.L + 1
        self.name = name
        if not name:
            try:
                self.name = tseries.name
            except AttributeError:
                pass
        
        # Standardize each time series
        self.standardize = True
        self.std_TS = self.standardized()
        
        # Embed the time series in a trajectory matrix: the first K columns
        # contain the first time series, the second K columns contain the
        # second time series, etc.
        self.nseries = self.std_TS.shape[1]
        self.X = np.array([self.std_TS.values[i:self.L+i, j] \
                           for j in range(0, self.nseries) \
                           for i in range(0, self.K)]).T

        # Decompose trajectory matrix using Singular Value Decomposition (SVD)
        self.U, self.Sigma, VT = np.linalg.svd(self.X)
        self.d = np.linalg.matrix_rank(self.X)

        # Reconstructed time series matrix (to fill)
        self.TS_comps = np.zeros((self.N*self.nseries, self.d))

        if not save_mem:
            # Construct and save all the elementary matrices
            self.X_elem = np.array([self.Sigma[i]*np.outer(self.U[:,i], VT[i,:])
                                    for i in range(self.d)])

            # Diagonally average the elementary matrices, store them as columns
            # in array.
            ind = np.arange(0, self.X.shape[1]+1, self.X.shape[1]//self.nseries)
            ts_temp = np.array_split(self.TS_comps, self.nseries, axis=0)
            for k in range(self.nseries):
                for i in range(self.d):
                    X_rev = self.X_elem[i, ::-1]
                    ts_temp[k][:,i] = \
                        [X_rev[:,ind[k]:ind[k+1]].diagonal(j).mean() \
                        for j in range(-X_rev.shape[0]+1, X_rev.shape[1]//2)]
            self.TS_comps = np.concatenate(ts_temp, axis=0)

            # Save the (transposed) V array
            self.V = VT.T

        else:
            # Reconstruct the elementary matrices without storing them
            ind = np.arange(0, self.X.shape[1]+1, self.X.shape[1]//self.nseries)
            ts_temp = np.array_split(self.TS_comps, self.nseries, axis=0)
            for k in range(self.nseries):
                for i in range(self.d):
                    X_elem = self.Sigma[i] * np.outer(self.U[:,i], VT[i,:])
                    X_rev = X_elem[::-1]
                    ts_temp[k][:,i] = \
                        [X_rev[:,ind[k]:ind[k+1]].diagonal(j).mean() \
                        for j in range(-X_rev.shape[0]+1, X_rev.shape[1]//2)]
            self.TS_comps = np.concatenate(ts_temp, axis=0)

            self.X_elem = 'Re-run with save_mem=False to retain the '\
                          'elementary  matrices.'

            # The V array may also be very large under these circumstances, so
            # discard it.
            self.V = 'Re-run with save_mem=False to retain the V matrix.'
        
        # Calculate the w-correlation matrix.
        self.calc_wcorr(w_repeats=self.nseries)

        # Split the time series and store in a new dimension and get the
        # dimensions to work with plotting methods.
        splits = np.array_split(self.TS_comps, self.nseries, axis=0)
        self.TS_comps = np.swapaxes(np.array(splits), 1,2).T

    @staticmethod
    def subplot_dims(k):
        """Return the numbers of rows and columns in which to arrange time 
        series subplots."""
        if k < 3:
            return k, 1
        else:
            return int(np.ceil(k/2)), 2

    def generate_subplots(self, k, row_wise=False):
        """Make separate subplots for each time series when plotting."""
        nrow, ncol = self.subplot_dims(k=self.nseries)
        fig, axes = plt.subplots(nrow, ncol, sharex=True)

        # Check if it's an array. If there's only one plot, it's just an Axes obj.
        if not isinstance(axes, np.ndarray):
            return fig, [axes]
        else:
            # Set the traversal: 'F' is col-wise, 'C' is row-wise
            axes = axes.flatten(order=('C' if row_wise else 'F'))
            # Delete any unused axes from the figure so that they don't show blank x- and y-axis lines
            for idx, ax in enumerate(axes[k:]):
                fig.delaxes(ax)
                # Turn ticks on for the last ax in each column, wherever it lands
                idx_to_turn_on_ticks = idx + k - ncol \
                                       if row_wise else idx + k - 1
                for tk in axes[idx_to_turn_on_ticks].get_xticklabels():
                    tk.set_visible(True)
            axes = axes[:k]
            return fig, axes
                                    
    def plot_all_patterns(self, legend=True):
        # Fiddle with color cycle - need more colors!
        plt.rcParams['font.size'] = 14
        plt.rcParams['image.cmap'] = 'plasma'
        plt.rcParams['axes.linewidth'] = 2
        color_cycle = cycler(color=plt.get_cmap('tab20').colors)
        
        fig, axes = self.generate_subplots(k=self.nseries)
        for i, ax in enumerate(axes):
            ax.plot(self.get_original().iloc[:,i], alpha=1, lw=1,
                             color='black', linestyle='dashed', zorder=0)
            for key in self.patterns.keys():
                ax.plot(self.reconstruct(self.patterns[key])\
                                 .iloc[:,i], lw=2, zorder=1)
            ax.set_xlabel("$t$")
            ax.set_ylabel(r"$\tilde{F}_i(t)$")
            if legend:
                legend = [r"$\tilde{F}_{%s}$" %(i-1)
                          for i in range(len(self.patterns.keys())+1)] + ["$F$"]
                ax.legend(legend, loc=(1.05,0.1))        
        fig.suptitle("{} Patterns of Drifter {}"\
                       .format(len(self.patterns), self.name))
        fig.tight_layout()
        plt.show()

    def plot_pattern(self, pattern, legend=True):
        # Fiddle with color cycle - need more colors!
        plt.rcParams['font.size'] = 14
        plt.rcParams['image.cmap'] = 'plasma'
        plt.rcParams['axes.linewidth'] = 2
        color_cycle = cycler(color=plt.get_cmap('tab20').colors)

        fig, axes = self.generate_subplots(k=self.nseries)
        for i, ax in enumerate(axes):
            ax.plot(self.get_original().iloc[:,i], alpha=1, color='black', lw=1,
                    linestyle='dashed', zorder=0, label='Original')
            ax.plot(self.reconstruct(self.get_patterns()[pattern]).iloc[:,i],
                    label='Pattern', linestyle='dotted')
            for v in self.get_patterns()[pattern]:
                ax.plot(self.reconstruct(int(v)).iloc[:,i],
                        label='Mode {}'.format(v), zorder=1)
            ax.set_xlabel("$t$")
            ax.set_ylabel(r"$\tilde{F}_i(t)$")
            if legend:
                ax.legend(loc=(1.05,0.1))
        fig.tight_layout()
        plt.show()
                    
    def __str__(self):
        return 'MSSA object containing multivariate singular spectrum '\
               'analysis results: {}'\
               .format(', '.join(list(self.__dict__.keys())))
    
    def __repr__(self):
        return '{lag_window:' + str(self.L) + ', tseries_length:' + \
                str(self.N) + ', trajMatLen:' + str(self.K)+'}'

def data_readIn(file):
    """
    Read in drifter data.
    """
    if not args.days:
        x = pd.concat([pd.read_csv(f, index_col=0) for f in files],
                       ignore_index=True)
    else:
        x = pd.concat([pd.read_csv(f, delimiter='\t', index_col='datetime',
                       parse_dates=True) for f in files], ignore_index=True)
    return x.dropna(how='all', axis=0)

def format_tseries(x, idx, splits):
    """
    Format time series for SSA: extract desired length and remove leading NAs.
    """
    ts_splits = np.array_split(x, splits)
    ts = pd.concat([i.reset_index(drop=True).dropna() for i in ts_splits],
                   axis=1)
    ts.columns = ['TS{}'.format(str(i)) for i in range(len(ts_splits))]
    ts.index = idx[:len(ts)]
    return ts
    
def normalize(x):
    """Normalize to interval [0,1]."""
    mins = x.min()
    maxs = x.max()
    return (x - mins) / (maxs - mins)
        
def ktest(tseries, m, k_thresh):
    """
    Finds values of k such that the first k elementary matrices from SSA 
    capture the desired variance given by 'k_thresh'. k corresponds to the 
    number of elementary components to sum to generate a reconstructed time 
    series.
    
    Returns a tuple containing k and the corresponding mean squared error 
    between the reconstructed time series and 'tseries'.

    Parameters
    ----------
    tseries : A list, 1D array, or pandas Series or DataFrame containing one or 
        more a time series on which to perform SSA. If more than one time 
        series is passed, a multivariate SSA will performed; otherwise, a
        standard univariate SSA.
    m : An integer representing the number of lags for creating time series 
        matrix.
    k_thresh : Float 0.0 < k_thresh < 1.0 specifying the fraction of total 
        variance to capture in the reconstructed time series.
    """
    # Perform SSA or MSSA
    if len(tseries.shape) == 1:
        ssa = SSA(tseries, window=m)
    else:
        ssa = MSSA(tseries, window=m)
    
    # Get k that satisfies k_thresh
    k = list(ssa.get_variance() > k_thresh).index(True)
    
    # Calculate mean squared reconstruction error
    err = ssa.calc_msre(slice(None, k))
    
    return (k, err)

def mtest(tseries, m_list, k_thresh):
    """
    Find optimal lag window m in m_list using ktest() with 'k_thresh'.
    
    Returns dictionary containing best m, k, and corresponding mean squared
    error between reconstructed time series and 'tseries'.
    
    Parameters
    ----------
    tseries : A list, 1D array, or pandas Series containing a time series for 
        SSA.
    m_list : A list containing the number of lags to test for creating time 
        series matrix.
    k_thresh : Float 0.0 < k_thresh < 1.0 specifying the fraction of total 
        variance to capture in the reconstructed time series.
    """
    # Check for time series length being less than expected.
    if sum(m_list > len(tseries)/2) > 0:
        try:
            print('Warning: Time series named {} is shorter than desired '\
                  'length. Some values of m will not be checked.'\
                  .format(tseries.name))
            print('Time series name: {}'.format(tseries.name))
        except AttributeError:
            print('Warning: Time series is shorter than desired length. '\
                  'Some values of m will not be checked.')
        # Use only m values that satisfy the condition 2 <= m <= len(tseries)/2
        m_list = m_list[m_list <= len(tseries)/2][:]
    
    results = {'TimeSeries': tseries}
    for m in m_list:
        k, err = ktest(tseries, m, k_thresh)
        if 'error' not in results:
            results['error'] = err
            results['k'] = k
            results['m'] = m
        elif err < results['error']:
            results['error'] = err
            results['k'] = k
            results['m'] = m
    return results

def best_ssa(tseries, m_list, k_thresh):
    """
    Creates an instance of class 'SSA' using the best lag window m and k 
    elementary components for time series reconstruction.
    
    Parameters
    ----------
    tseries : A list, 1D array, or pandas Series containing a time series for 
        SSA.
    m_list : A list containing the number of lags to test for creating time 
        series matrix.
    k_thresh : Float 0.0 < k_thresh < 1.0 specifying the fraction of total 
        variance to capture in the reconstructed time series.
    """
    m_results = mtest(tseries=tseries, m_list=m_list, k_thresh=k_thresh)
    return SSA(tseries=tseries, window=m_results['m'])

def signal_noise_split(m, method=['counts', 'pdf'], cor_thresh=1.0e-2, var_thresh=3.0e-3, print_var=False):
    """
    Returns the splitting row/column index (of the correlation matrix) between
    signal and noise. Two methods are supported:
    
    'counts' determines for each row the number of correlations greater than 
    the correlation threshold. Those rows containing greater than 30% of the 
    total number of elements (rows/columns) will be considered noise.
    
    'pdf' uses the variance of the PDF accross rows/columns of the
    correlation matrix. A PDF is calculated on the correlation matrix (ignoring 
    the diagonal of perfect correlation) and the variance for each row is 
    calculated. This list of variances is checked for the FIRST element that 
    exceeds the variance threshold and the index of this element is returned.
    All rows/columns beyond the first row/column that exceeds the variance
    threshold is taken as noise.
    
    Parameters
    ----------
    m : A correlation matrix or array.
    method: A string, either 'counts' or 'pdf', indicating the splitting method 
        to use. Defaults to 'counts'.
    cor_thresh : A float indicating the maximum correlation to be considered 
        significant, used with method='counts'. All values less than this will
        be assumed no correlation. Defaults to 0.01.
    var_thresh : A float indicating the variance threshold to split the
        correlation matrix, used with method='pdf'. All rows/columns with PDF
        variances less than this threshold should be taken as signal, while all
        rows/columns with variances greater than this threshold should be taken
        as noise. Defaults to 0.003.
    print_var : Logical indicating whether variance of the PDF function should 
        be printed, usually for diagnostic purposes. Only works for 
        method='pdf'. Defaults to False.
    """
    assert method in ['pdf', 'counts', ['counts','pdf'], ['pdf', 'counts']],\
        "Invalid method option. Must be either 'counts' or 'pdf'."
    n = m.shape[0]
    # Set default method of both options are passed
    if isinstance(method, list):
        method = 'counts'
    if method.casefold() == 'pdf':
        # Use NumPy as_strided to push diagonal to the first column, slice it
        # off, and reshape to proper number of columns. This technique is
        # adapted from https://stackoverflow.com/questions/29394377/minimum-of-
        #    numpy-array-ignoring-diagonal
        no_diag = as_strided(m, (n-1,n+1), (m.itemsize*(n+1),
                             m.itemsize))[:,1:].reshape(n,n-1)
        # Calculate variances across columns
        var = norm.cdf(no_diag).std(1)**2
        if print_var:
            print(var)
        # Return first index where variance exceeds supplied threshold
        return next(np.where(var == i) for i in var if i >= var_thresh)[0][0]
    else:
        def count_first_True(a):
            """Count the first consecutive True bools, stopping at the first 
            False and ignoring any subsequent True bools."""
            maxidx = a.argmax()
            pos = a[maxidx:].argmin()
            if a[maxidx]:
                if pos==0:
                    return a.size - maxidx
                else:
                    return pos
            else:
                return 0
        counts = (m>=cor_thresh).sum(0) < 0.3*n
        return count_first_True(counts)

def perform_mssa(m=None, rm_mean=False, norm=False, silent=False):
    # Make time series
    ts = format_tseries(df.loc[:,drifter], idx=index.ptime, splits=len(files))

    # Subtract mean (optional)
    if rm_mean:
        ts = ts - ts.mean()
    
    # Normalize (optional)
    if args.normalize:
        ts = normalize(ts)
    
    # Decompose
    df_ssa = MSSA(tseries=ts, window=m)
    if not silent:
        print('Performing SSA on drifter ID {} with window of size {} lags '\
              '(of max {}).'.format(drifter, m, len(ts)//2))
    return df_ssa


# =========================================================================== #
# Read in data
args = parse_args()
files = list(args.file)
df = data_readIn(files)
index = pd.read_csv(files[0], usecols=[0])

# Empty dataframes to store results
df_summary = pd.DataFrame(index=df.columns, columns=['LagWindow', 'TSeriesLen',
                         'TrajMatLen', 'SignalComponents', 'Patterns'])
df_summary.index.name = 'drifter'

# Perform MSSA
if not args.days:
    m_wind = int(args.lag * df.shape[0]//len(files))
else:
    m_wind = int(args.lag * 4*24 * args.days)

for drifter in df.columns:
    df_ssa = perform_mssa(m=m_wind, norm=args.normalize, silent=True)
    df_summary.loc[drifter] = df_ssa.find_patterns(signal_cutoff=args.signal,
                                                 corr_strength=args.correlation)

# Write to file or print summary
if args.writeout:
    prefix = os.path.splitext(os.path.basename(files[0]))[0].split('_')[0]
    print('Pattern counts:\n')
    print(df_summary.Patterns.value_counts())
    print('Written to {}_mssa_patterns.csv'.format(prefix))
    df_summary.to_csv('{}_mssa_patterns.csv'.format(prefix))
else:
    print('Pattern counts:\n')
    print(df_summary.Patterns.value_counts())

# =========================================================================== #
