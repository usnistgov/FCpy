# -*- coding: utf-8 -*-
"""
NIST-developed software is provided by NIST as a public service. You may use, 
copy, and distribute copies of the software in any medium, provided that you 
keep intact this entire notice. You may improve, modify, and create derivative 
works of the software or any portion of the software, and you may copy and 
distribute such modifications or works. Modified works should carry a notice 
stating that you changed the software and should note the date and nature of 
any such change. Please explicitly acknowledge the National Institute of 
Standards and Technology as the source of the software. 

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY 
OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, 
INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER 
REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED 
OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR 
MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS 
THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, 
OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and 
distributing the software and you assume all risks associated with its use, 
including but not limited to the risks and costs of program errors, compliance 
with applicable laws, damage to or loss of data, programs or equipment, and the 
unavailability or interruption of operation. This software is not intended to 
be used in any situation where a failure could cause risk of injury or damage 
to property. The software developed by NIST employees is not subject to 
copyright protection within the United States.


@author: Evan Groopman
@affiliation: National Institute of Standards and Technology (NIST)
@email: evan.groopman@nist.gov

Generate Feldman-Cousins (1998) tables and interpolate to create CLs and ULs
https://journals.aps.org/prd/abstract/10.1103/PhysRevD.57.3873

Pure Python & numpy implementation.


For Poisson statistics with background
Known measured background and measured events.

Feldman and Cousins generate their tables by calculating mu on the interval [0,50]
in steps of 0.005, and b in the interval [0,25] in steps of 0.001 (and some additional 
searching to finer precision). n being discrete results in mild pathologies,
see paper.

Add correction for low counts put forth by Roe 1999.
Replaces Poisson probability density with conditional probability

Benchmarks:

"""

import numpy as np
from scipy.stats import norm
from scipy.special import gammaln #factorial
import warnings

try:
    import dask
except ModuleNotFoundError:
    dask = None


#FeldmanCousins(1998) Poisson definitions
def poissonPMF(k, mu):
    """
    Poisson probability mass function (pmf). This was the fastest implementation.

    Parameters
    ----------
    k : int, np.ndarray
        Number of Poisson events
    mu : float, np.ndarray
        Mean Poisson rate

    Returns
    -------
    Vectorized probability for each k,mu.

    """
    # 
    #using gammaln is much faster than factorial function
    #need to make sure log and exp return the proper values when given mu is 0
    return(np.exp(k * np.log(mu, out=np.zeros_like(mu, dtype=float), where=mu!=0.0) - mu - gammaln(k + 1), 
                  out=np.exp(-mu), where=mu!=0.0))

def Pn_mu(mu, b, n):
    """Wrapper for poissonPMF(k, mu) to include known background rate, b.
    

    Parameters
    ----------
    mu : float, np.ndarray
        Mean Poisson rate
    b : float, np.ndarray
        Known background rate
    n : int, np.ndarray
        Number of Poisson events

    Returns
    -------
    Poission probability at n events with rate mu+b.

    """
    return(poissonPMF(n, mu+b))

def R_poisson(mu, b, n):
    """
    Likelihood ratio P(n,mu)/P(n,mu_best).
    Feldman and Cousins (1998) use this ratio for their ordering principle

    Parameters
    ----------
    mu : float, np.ndarray
        True count rate.
    b : float, np.ndarray
        Mean background rate.
    n : int, np.ndarray
        Number of observed counts.

    Returns
    -------
    Likelihood ratio for ordering.

    """
    return(Pn_mu(mu, b, n)/Pn_mu(np.maximum(0, n - b), b, n))

#corrections based upon Roe & Woodroofe (1999). Poisson case only
def _qn_mu_lt(mu, b, n, n0):
    """internal function for k < n0"""
    return(poissonPMF(n, mu+b)/np.sum(poissonPMF(np.arange(0,n0+1), np.ones(n0+1, dtype=float)*b)))

def _qn_mu_gt(mu, b, n, n0):
    """internal function for k > n0"""
    #need to add a third dimension and then collapse back to 2
    #needs a sum over j=0 to n0 of p_b(j)*p_mu(k-j)
    j = np.arange(0,n0+1)
    denom = np.sum(poissonPMF(j, b*np.ones(n0+1, dtype=float)))
    
    #b is a scalar, n0 is a scalar
    #mu is a 2d array, n is a 2d array from meshgrid
    # MU = np.repeat(mu[np.newaxis,:,:], n0+1, axis=0) #repeat mu into z direction
    nn = np.broadcast_to(n, (n0+1, *n.shape)) - j[:,np.newaxis,np.newaxis] #k-j in Roe1999
    mm = np.broadcast_to(mu, (n0+1, *mu.shape))
    
    num = (poissonPMF(j, b*np.ones(n0+1, dtype=float))[:,np.newaxis,np.newaxis] * poissonPMF(nn, mm)).sum(axis=0)
    
    return(num/denom)

def qn_mu(mu, b, n, n0):
    """calculates Roe & Woodroofe (1999) qn_mu correction."""
    return(np.select([n <= n0, n > n0], [_qn_mu_lt(mu,b,n,n0), _qn_mu_gt(mu,b,n,n0)]))

def Rq_poisson(mu, b, n, n0):
    """Likelihood ratio Q(n,mu)/Q(n,mu_best).
    Roe & Woodroofe use this ratio for their ordering principle"""
    return(qn_mu(mu, b, n, n0)/qn_mu(np.maximum(0, n - b), b, n, n0))

#Feldman-Cousins(1998) gaussian. #needs a tiny bit of work and optimization
def gaussian(x, mu):
    """Gaussian pdf"""
    return(np.exp(-(x - mu)**2/2)/np.sqrt(2*np.pi))

def _R_gauss_pos(x, mu):
    return(np.exp(-(x - mu)**2/2))

def _R_gauss_neg(x, mu):
    return(np.exp(x*mu - mu**2/2))

def R_gauss(x, mu):
    """when x>=0 use R_gauss_pos(); when x<0 use R_guass_neg().
    Gaussian with cutoff"""
    return(np.select([x >= 0, x < 0], [_R_gauss_pos(x, mu), _R_gauss_neg(x, mu)]))


def FC_gauss(x0, conf=0.95):
    """Feldman Cousins confidence intervals for the Gaussian case with a lower bound.
    
    only use for small x0 < 10!"""
    #this is not optimized, but it works and relatively quickly. FC_poisson is much faster with better searching
    #based on the definition of gaussian() above, sigma = 1
    sigma = 1
    
    xmax = np.clip(x0 + 5*sigma, 3, 10)
    xmin = np.clip(x0 - 4*sigma, -10, None)
    # xmax = 5
    # xmin = -4
    roughstep = 0.01
    nsteps = np.maximum(10, int((xmax-xmin)/roughstep))
    xx = np.linspace(xmin, xmax, nsteps)
    #add in x0 so it's always there in xx
    if x0 not in xx:
        xx = np.hstack([xx,x0])
        xx.sort()
    xind = np.where(xx==x0)[0].item()
    
    mumin = 0
    mumax = x0 + 4*sigma #abs(xmax) #abs(xmin) + abs(xmax) #same range as x, but strictly positive
    murange = np.linspace(mumin, mumax, nsteps)
    
    XX, MU = np.meshgrid(xx, murange)
    
    pdf = gaussian(XX, MU)
    # #calculate R at each point
    r = R_gauss(XX, MU)
    ind = np.argsort(r, axis=1)[:,::-1]

    cutoff_ind = np.sum(np.cumsum(pdf[np.repeat(np.arange(len(murange))[:,np.newaxis], len(xx), axis=1), ind], axis=1) <= conf/roughstep, axis=1)
    ranks = np.argsort(ind, axis=1) #much faster than scipy.stats.rankdata
    grid = ranks <= cutoff_ind[:, np.newaxis]
    
    #if the whole line is true, this can mess up the next min/max step
    grid[np.all(grid, axis=1)] = False
    #if the last index in x is true, can mess up the min/max
    grid[grid[:,-1]] = False
    gridline = grid[:,xind]
    
    muline = murange[gridline]
    CL_low_rough = muline.min()
    CL_high_rough = muline.max()

    return(CL_low_rough, CL_high_rough)

   

def FC_poisson(n0, b, t, conf=0.95, useCorrection= False, tol=5E-4,
               mumin= None, mumax= None,
               fixPathology= False, bRange=0.01, bStep= 0.001):
    """
    Feldman Cousins confidence level calculation for the Poisson case with known mean background.
    
    Uses numpy meshgrid to get all permutations of n and mu(+b) to avoid explicit loops.
    Much much faster than for/while loops in Python.

    Parameters
    ----------
    n0 : int
        Number of observed counts
    b : float
        Mean background count rate (counts/s).
    t : float
        Total measurement time (s).
    conf : float between [0,1.0], optional
        Confidence level to calulate. The default is 0.95.
    useCorrection : bool, optional
        Use the Roe & Woodroofe (1999) correction for low counts. The default is False.
        This is computationally more epxensive, especially for n0 > ~30
    tol : float, optional
        Calculation tolerance for mu. The default is 5E-4.
    mumin : float or None, optional
        If given, provide the minimum mu value to search. Otherwise, one selected automatically.
    mumax : float or None, optional
        If given, provide the maximum mu value to search. Otherwise, one selected automatically.
    fixPathology : bool, optional
        Feldman-Cousins noted that a mild pathology occurs for fixed n0 with varying b.
        This searches the neighborhood of b and ensures that the upper CI limit is 
        monotonically decreasing as a function of b across the search range.
    bRange : float, optional
        Range over which to search for pathology np.clip([b - bRange, b + bRange], 0, None)
    bStep : float, optional
        Step size to use in pathology search

        
    Returns
    -------
    np.ndarray([CI lower limit, CI upper limit])

    """
    
    if fixPathology:
        
        bRange = abs(bRange); bStep = abs(bStep)
        if bRange < bStep:
            bRange = bStep
            
        # bMin, bMax = np.clip([b - bRange, b + bRange], 0, None)
        # nSteps = int((bMax-bMin)/bStep)
        # Bs = np.linspace(bMin, bMax, num=nSteps, endpoint=True)
        # insert real value at index
        # ind = np.searchsorted(Bs, b)
        # Bs = np.insert(Bs, ind, b)
        nSteps = np.maximum(int(bRange/bStep), 2) #ensure at least 2 steps
        Bs = np.linspace(b, b + bRange, num=nSteps, endpoint=True)
        
        if dask is not None:
            results = []
            for B in Bs:
                res = dask.delayed(_FC_poisson)(n0, B, t, conf, useCorrection, tol, mumin, mumax, False)
                results.append(res)
            results = np.array(dask.compute(*results))
            
        else:
            results = np.empty((nSteps, 2))
            for i, B in enumerate(Bs):
                results[i,:] = _FC_poisson(n0, B, t, conf=conf, useCorrection=useCorrection, tol=tol,
                                           mumin= mumin, mumax= mumax,
                                           useDask= False)
        
        CI = results[0,:]
        
        #check if monotonically decreasing after b
        diffCheck = (results[1:, 1] - results[:-1,1]) > 0
        if np.any(diffCheck):
            ind = np.where(diffCheck)[0]
            CIup = results[ind][1]
            if CIup > CI[1]:
                CI[1] = CIup
        return(CI)
        
    else:
        #swich on dask internally if high enough counts or RW correction is on.
        if (n0 + b > 3) or useCorrection:
            useDask = True
        else:
            useDask = False
        return(_FC_poisson(n0, b, t, conf=conf, useCorrection=useCorrection, tol=tol, 
                           mumin= mumin, mumax= mumax, useDask= useDask))
    
        
def _FC_poisson(n0, b, t, conf=0.95, useCorrection= False, tol=5E-4,
                 mumax= None, mumin= None, useDask= True):
    """
    Feldman Cousins confidence level calculation for the Poisson case with known mean background.
    
    Uses numpy meshgrid to get all permutations of n and mu(+b) to avoid explicit loops.
    Much much faster than for/while loops in Python.

    Parameters
    ----------
    n0 : int
        Number of observed counts
    b : float
        Mean background count rate (counts/s).
    t : float
        Total measurement time (s).
    conf : float between [0,1.0], optional
        Confidence level to calulate. The default is 0.95.
    useCorrection : bool, optional
        Use the Roe & Woodroofe (1999) correction for low counts. The default is False.
        This is computationally more epxensive, especially for n0 > ~30
    tol : float, optional
        Calculation tolerance for mu. The default is 5E-4.
    mumin : float or None, optional
        If given, provide the minimum mu value to search. Otherwise, one selected automatically.
    mumax : float or None, optional
        If given, provide the maximum mu value to search. Otherwise, one selected automatically.
    useDask : bool, optional
        Use dask package, if available, to speed up some parallel computations
        
    Returns
    -------
    np.ndarray([CI lower limit, CI upper limit])

    """

    n0 = int(n0)
    
    bt = b*t
    
    sigma = abs(norm.ppf((1-conf)/2))
    

    
    # mumax = n0+b + 5*np.sqrt(n0+b)
    if n0 + bt <= 3:
        if mumax is None:
            mumax = 11.0
        if mumin is None:
            mumin = 0.0
    else:
        #scale by the desired conf?
        if mumax is None:
            mumax = n0 + bt + 5*np.sqrt(n0+bt)
        if mumin is None:
            mumin = np.maximum(0, n0 - 5*np.sqrt(n0)) #when looking for mu lower limit, need to ignore contribution from b
    
    #first do a coarse scan to find approximate locations of CLs
    #then do smaller scans near limits with precision of tol
    #minimum 10 initial steps
    roughstep = 200*tol #0.1 #this seems to work well
    steps = np.maximum(10, int((mumax - mumin)/roughstep))
    
    murange = np.linspace(mumin, mumax, steps)
    
    if n0 <= 100:
        minn = 0
        maxn = (2*n0 + 16 + bt)
        if conf > 0.95:
            maxn *= (sigma/1.96) #is later rounded to nearest int
            maxn = int(maxn)
    else:
        minn = 0
        #this seems to work
        maxn = 1.75*n0 + bt
    
    nn = np.arange(minn, maxn, dtype=int) #could collapse further using estimated minn
    
    #if n is large enough, can really collapse this grid by using a Gaussian approximation of the CLs
    #CL_low_rough will be ~ n0 - sqrt(n0)*sigma + 1
    #CL_high_rough will be ~ n0 + sqrt(n0)*sigma + 1
    #remove central murage that is pretty far outside this area. no need to calculate it
    if n0 >= 50: #set 50 as the limit for switching to Gaussian approx and reducing rough scan area
        CL_low_rough = n0 - np.sqrt(n0)*sigma + 1 + n0*roughstep #set bar above CL_rough and remove in between rough estimates
        CL_high_rough = n0 + np.sqrt(n0)*sigma + 1 - n0*roughstep
        murange = murange[(murange < CL_low_rough) | (murange > CL_high_rough)]
    
    
    NN, MU = np.meshgrid(nn, murange)
    
    if useCorrection:
        if bt > 0:
            pdf = qn_mu(MU, bt, NN, n0)
            r = Rq_poisson(MU, bt, NN, n0)
        else:
            pdf = qn_mu(MU, 1E-7, NN, n0)
            r = Rq_poisson(MU, 1E-7, NN, n0)

    else:
        pdf = Pn_mu(MU, bt, NN)
        # #calculate R at each point
        r = R_poisson(MU, bt, NN)
    
    ind = np.argsort(r, axis=1)[:,::-1]

    cutoff_ind = np.sum(np.cumsum(pdf[np.repeat(np.arange(len(murange))[:,np.newaxis], len(nn), axis=1), ind], axis=1) <= conf, axis=1)
    ranks = np.argsort(ind, axis=1) #much faster than scipy.stats.rankdata
    grid = ranks <= cutoff_ind[:, np.newaxis]
    
    
    #if the whole line is true, this can mess up the next min/max step
    grid[np.all(grid, axis=1)] = False
    gridline = grid[:,n0]
    
    muline = murange[gridline]
    CL_low_rough = muline.min()
    CL_high_rough = muline.max()
    # print(CL_low_rough, CL_high_rough)
    
    
    if dask is not None and useDask:
        CL_low = dask.delayed(_getLowerCI)(CL_low_rough, n0, bt, nn, conf, 
                              roughstep, tol, useCorrection)
        CL_high = dask.delayed(_getUpperCI)(CL_high_rough, n0, bt, nn, conf, 
                              roughstep, tol, useCorrection)
        results = dask.compute(CL_low, CL_high)
        
        CL_low, CL_high = results
    else:
        
        CL_low = _getLowerCI(CL_low_rough, n0, bt, nn, conf= conf, 
                              roughstep= roughstep, tol= tol, 
                              useCorrection= useCorrection)
        CL_high = _getUpperCI(CL_high_rough, n0, bt, nn, conf= conf, 
                              roughstep= roughstep, tol= tol, 
                              useCorrection= useCorrection)
    
    #raise warning if CL_low and CL_high are at the limits of mumax and mumin
    if CL_high >= mumax:
        warnings.warn('The upper limit is constrained by mumax = {} and may not be correct. Increase mumax or set to None.'.format(mumax),
                      category=RuntimeWarning)
    if CL_low <= mumin and mumin > 0:
        warnings.warn('The lower limit is constrained by mumin = {} and may not be correct. Lower mumin or set to None.'.format(mumin),
                      category=RuntimeWarning)

    return(np.array([CL_low, CL_high]))

def _getLowerCI(CL_low_rough, n0, bt, nn, conf= 0.95, roughstep= 0.1, tol= 5E-4, useCorrection= False):
    """Internal function for FC_poisson"""
    if n0 == 0 and bt == 0:
        CL_low = 0.0
    else:
        #do a fine loop to find the lower limit (steps of tol).
        mumin = np.max(np.array([0, CL_low_rough - 1.5*roughstep]))
        mumax = CL_low_rough + 1.5*roughstep
        murange = np.linspace(mumin, mumax, int((mumax - mumin)/tol))
        NN, MU = np.meshgrid(nn, murange)
        
        if useCorrection:
            if bt > 0:
                pdf = qn_mu(MU, bt, NN, n0)
                r = Rq_poisson(MU, bt, NN, n0)
            else:
                pdf = qn_mu(MU, 1E-7, NN, n0)
                r = Rq_poisson(MU, 1E-7, NN, n0)
            
        else:
            pdf = Pn_mu(MU, bt, NN)
            # #calculate R at each point
            r = R_poisson(MU, bt, NN)
        ind = np.argsort(r, axis=1)[:,::-1]
        
        cutoff_ind = np.sum(np.cumsum(pdf[np.repeat(np.arange(len(murange))[:,np.newaxis], len(nn), axis=1), ind], axis=1) <= conf, axis=1)
        ranks = np.argsort(ind, axis=1) #much faster than scipy.stats.rankdata
        grid = ranks <= cutoff_ind[:, np.newaxis]
        
        CL_low = np.min(murange[grid[:,n0]])
        if abs(1 - CL_low/tol) <= 0.005:
            CL_low = 0.0
        
    return(CL_low)

def _getUpperCI(CL_high_rough, n0, bt, nn, conf= 0.95, roughstep= 0.1, tol= 5E-4, useCorrection= False):
    """Internal function for FC_poisson"""
    #do a fine loop to find upper limit (steps of tol)
    mumin = np.max(np.array([0, CL_high_rough - 1.5*roughstep]))
    mumax = CL_high_rough + 1.5*roughstep
    murange = np.linspace(mumin, mumax, int((mumax - mumin)/tol))
    NN, MU = np.meshgrid(nn, murange)
    
    if useCorrection:
        if bt > 0:
            pdf = qn_mu(MU, bt, NN, n0)
            r = Rq_poisson(MU, bt, NN, n0)
        else:
            pdf = qn_mu(MU, 1E-7, NN, n0)
            r = Rq_poisson(MU, 1E-7, NN, n0)
    else:
        pdf = Pn_mu(MU, bt, NN)
        # #calculate R at each point
        r = R_poisson(MU, bt, NN)
    ind = np.argsort(r, axis=1)[:,::-1]
    
    cutoff_ind = np.sum(np.cumsum(pdf[np.repeat(np.arange(len(murange))[:,np.newaxis], len(nn), axis=1), ind], axis=1) <= conf, axis=1)
    ranks = np.argsort(ind, axis=1) #much faster than scipy.stats.rankdata
    grid = ranks <= cutoff_ind[:, np.newaxis]
    
    CL_high = np.max(murange[grid[:,n0]])
    return(CL_high)

def FC_poisson_confbands(nrange= [0,10], b= 0.0, t= 1, conf= 0.95, useCorrection= False, tol= 5E-4,
                         useDask= True):
    """
    Generate upper and lower confidence bands for the range of given n values

    Parameters
    ----------
    nrange : tuple, list, np.ndarray of ints
        Range over which to generate confidence bands. The default is [0,10].
    b : float, optional
        Average background rate. The default is 0.0.
    t : int, optional
        Measurement time in seconds. Used to generate the average number of background counts. The default is 1.
    conf : float, optional
        Confidence interval in the range [0,1]. The default is 0.95.
    useCorrection : bool, optional
        Use Roe&Woodroofe(1999) correction. The default is False.
    tol : float, optional
        Numerical tolerance on the confidence band. The default is 5E-4.

    Returns
    -------
    Nx2 np.ndarray where N is the number of integers in nrange

    """
    nrange = np.array(nrange).astype(int)
    nn = np.arange(nrange[0], nrange[1]+1, dtype=int)
    
    if dask is not None and useDask:
        results = []
        for n0 in nn:
            res = dask.delayed(_FC_poisson)(n0, b, t, conf, useCorrection, tol, False)
            results.append(res)
        results = np.array(dask.compute(*results))
        results = np.vstack([nn, results.T]).T
    else:
        results = np.empty((len(nn),3))
        for i, n0 in enumerate(nn):
            results[i,0] = n0
            results[i,1:] = _FC_poisson(n0, b, t, conf=conf, useCorrection=useCorrection, tol=tol, useDask= False)
    
    return(results)
    

def FC_poisson123(n0, b, t, useCorrection= False, tol=5E-4):
    """
    Convenience function for calculating 1-, 2-, and 3- sigma levels (Gaussian equivalent).
    conf = [0.683, 0.95, 0.997]. 

    Parameters
    ----------
    n0 : int
        Number of observed counts
    b : float
        Mean background rate.
    t : float
        Total measurement time.
    conf : float < 1.0, optional
        Confidence level to calulate. The default is 0.95.
    useCorrection : bool, optional
        Use the Roe & Woodroofe (1999) correction for low counts. The default is False.
        This is computationally more epxensive, especially for n0 > ~30
    tol : float, optional
        Calculation tolerance for mu. The default is 5E-4.

    Returns
    -------
    Tuple of CI arrays for conf = [0.683, 0.95, 0.997].

    """
    # 4 sigma = 0.9999 for reference
    
    n0 = int(n0)
    
    bt = b*t
    
    sigma = 3
    
    # mumax = n0+b + 5*np.sqrt(n0+b)
    if n0 + bt <= 3:
        mumax = 11.0
        mumin = 0.0
    else:
        #scale by the desired conf?
        mumax = n0 + bt + 5*np.sqrt(n0+bt)
        mumin = np.maximum(0, n0 - 5*np.sqrt(n0)) #when looking for mu lower limit, need to ignore contribution from b
    
    #first do a coarse scan to find approximate locations of CLs
    #then do smaller scans near limits with precision of tol
    #minimum 10 initial steps
    roughstep = 200*tol #0.1 #this seems to work well
    steps = np.maximum(10, int((mumax - mumin)/roughstep))
    
    murange = np.linspace(mumin, mumax, steps)
    
    if n0 <= 100:
        maxn = int((2*n0 + 16 + bt)*(sigma/1.96))
    else:
        maxn = 1.75*n0 + bt
    
    nn = np.arange(0, maxn, dtype=int)

    NN, MU = np.meshgrid(nn, murange)
    
    if useCorrection:
        if bt > 0:
            pdf = qn_mu(MU, bt, NN, n0)
            r = Rq_poisson(MU, bt, NN, n0)
        else:
            pdf = qn_mu(MU, 1E-7, NN, n0)
            r = Rq_poisson(MU, 1E-7, NN, n0)

    else:
        pdf = Pn_mu(MU, bt, NN)
        # #calculate R at each point
        r = R_poisson(MU, bt, NN)
    
    ind = np.argsort(r, axis=1)[:,::-1]

    ii = np.repeat(np.arange(len(murange))[:,np.newaxis], len(nn), axis=1)
    
    CIs= []
    ranks = np.argsort(ind, axis=1) #much faster than scipy.stats.rankdata
    for conf in [0.683, 0.95, 0.997]:
        cutoff_ind = np.sum(np.cumsum(pdf[ii, ind], axis=1) <= conf, axis=1)
        grid = ranks <= cutoff_ind[:, np.newaxis]
        #if the whole line is true, this can mess up the next min/max step
        grid[np.all(grid, axis=1)] = False
        gridline = grid[:,n0]

        muline = murange[gridline]
        CL_low_rough = muline.min()
        CL_high_rough = muline.max()

        CL_low = _getLowerCI(CL_low_rough, n0, bt, nn, conf= conf, 
                              roughstep= roughstep, tol= tol, 
                              useCorrection= useCorrection)
        CL_high = _getUpperCI(CL_high_rough, n0, bt, nn, conf= conf, 
                              roughstep= roughstep, tol= tol, 
                              useCorrection= useCorrection)
        
        CIs.append((CL_low, CL_high))
    
    return(CIs)

def FC_poisson_list(n= [0], b= 0.0, t= [1], conf= 0.95, useCorrection= False, tol= 5E-4,
                         useDask= True):
    """
    Generate CIs for a list of n, t, and conf values.

    Parameters
    ----------
    n : int, tuple, list, np.ndarray of ints
        Collection of n values to calculate. The default is [0].
    b : float, (optional)
        Average background rate. The default is 0.0.
    t : int, float, tuple, list, ndarray (optional)
        Measurement time in seconds. If an int/float, the same value is used for all n.
        Otherwise must be the same length as n.
        Used to generate the average number of background counts. The default is 1.
    conf : float, (optional)
        Confidence interval in the range [0,1]. The default is 0.95.
    useCorrection : bool, (optional)
        Use Roe&Woodroofe(1999) correction. The default is False.
    tol : float, (optional)
        Numerical tolerance on the confidence band. The default is 5E-4.

    Returns
    -------
    Nx2 np.ndarray where N is the number of n

    """
    
    if isinstance(n, (int, float)):
        n = int(n)
        if isinstance(t, (list, tuple, np.ndarray)):
            raise ValueError('cannot compute with multiple t values for a single n')
    elif isinstance(n, (list, tuple, np.ndarray)):
        n = np.asarray(n, dtype=int)
        if isinstance(t, (list, tuple, np.ndarray)):
            t = np.asarray(t, dtype=float)
            if len(t) != len(n):
                raise ValueError('n and t have different shapes: {} {}'.format(n.shape, t.shape))
        else:
            t = np.ones_like(n, dtype=float) * t

    if dask is not None and useDask:
        results = []
        for n0, t0 in zip(n,t):
            res = dask.delayed(_FC_poisson)(n0, b, t0, conf, useCorrection, tol)
            results.append(res)
        results = np.array(dask.compute(*results))
    else:
        results = np.empty((len(n),2))
        for i, (n0,t0) in enumerate(zip(n,t)):
            results[i,:] = _FC_poisson(n0, b, t0, conf=conf, useCorrection=useCorrection, tol=tol)
    
    return(results)
    
if __name__ == "__main__":
    #command line usage
    import sys
    
    # example: python -m FC 10 0 100 0.95 # OR
    # python FC.py 10 0 100 0.95
    
    args = sys.argv #n0, b, t, conf
    # print(sys.argv)
    fname = args[0]
    n0 = int(args[1])
    b = float(args[2])
    t = float(sys.argv[3])
    conf = float(sys.argv[4])
    
    LL, UL = FC_poisson(n0, b, t, conf=conf)
    sys.stdout.write(str('%.3f %.3f' % (LL, UL)))
    sys.exit(0)
