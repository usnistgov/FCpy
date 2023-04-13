# FCpy

Calculate Feldman-Cousins confidence intervals with Python & numpy. 

Highly vectorized and fast implementation of Poisson and Gaussian confidence interval (CI) calculation algorithms described in Feldman & Cousins (1998) Phys. Rev. D, 57, 3873 https://doi.org/10.1103/PhysRevD.57.3873. Generally for Poisson processes in the presence of a known mean background count rate or Gaussian processes with a physical lower cutoff of 0. Note, these do not incorporate uncertainty on the mean background rate (see Bayesian example instead).

Includes optional low-count correction to Poisson CIs described by Roe & Woodroofe (1999) Phys. Rev. D, 60, 053009 https://doi.org/10.1103/PhysRevD.60.053009 (this is more computationally expensive, see Benchmarks). The original Feldman & Cousins formulation would result in tighter confidence intervals with increasing background levels, which does not make sense. This is especially problematic for the case of 0 total counts, where the number of background counts were clearly 0. If, for example, the average background rate were relatively large, then this measurement would not be representative of the true mean background - it would be a statistically unlikley, yet possible event. Therefore, the confidence on the non-background counts should not be reduced - this measurement was clearly an outlier.

Alternatively, a Bayesian example shows the use of PyMC to perform the analogous MCMC simulation. The benefits of the Bayesian solution are that it can incorporate uncertainties on the background rate and other instrumental sensitivity factors. The Bayesian solution also yields the CI of the parameter of interest (e.g., the percent probability that the number of counts lies within the CI), which is generally what the user actually wants. The frequentist interpretation is instead that upon repeating the same experiment X times, the result will fall within that range CI % of the repetitions.

Further reading: see Coakley et al. (2010) Meas. Sci. Tech. 21, 035102 https://doi.org/10.1088/0957-0233/21/3/035102 for a discussion of coverage properties between FC, Bayesian, and propagation of errors (POE) methods.


## Usage:

```python
from FCpy.FCpy import FC
```

Function documentation:
```python
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
```

Calculate 95% Feldman-Cousins CI for 5 counts in 1000 s with mean background rate of 0.01 counts/s (no correction).
```python
FC.FC_poisson(5, 0.01, 1000, conf= 0.95, useCorrection= False)
# [Out]: array([0. , 2.93706109])
```

95% CI with Roe & Woodroofe correction.
```python
FC.FC_poisson(5, 0.01, 1000, conf= 0.95, useCorrection= True)
# [Out]: array([0. , 4.76047702])
```

A convenience function ```FC_poisson123()``` calculates 68.3%, 95%, and 99.7% CIs:

```python
FC.FC_poisson123(5, 0.01, 1000, useCorrection= False)
# [Out]: [(0.0, 0.27183253885023984), (0.0, 2.937061090328606), (0.0, 6.40224918079424)]
```

A convenience function ```FC_poisson_confbands()``` calculates confidence bands for a range of n (optionally uses dask to parallelize operations):

```python
FC_poisson_confbands(nrange= [0,10], b= 0.0, t= 1, conf= 0.95, useCorrection= False, tol= 5E-4,
                     useDask= True)
```

Command-line usage (defaults to Poisson CI without correction):

```shell
$ python -m FC 5 0.01 1000 0.95
```

OR

```shell
$ python FC.py 5 0.01 1000 0.95
```

There is also a small GUI written in Python and Qt for calculating Poisson confidence intervals.

## Confidence Bands

The image below shows a comparison of confidence bands using the Feldman-Cousins (FC), Roe-Woodroofe (RW), and Bayesian methods for different average detector background rates over a 400 s measurement. On the Large-Geometry Secondary Ion Mass Spectrometer at NIST, the long-term average electron multiplier detector background is 0.0013 counts/s. All three methods are in general agreement until the number of background counts is comparable to or larger than the signal of interest. Oddly, there is a pathology in the RW correction for n = 4 at this background level, but otherwise it finds good agreement with the Bayesian model. The FC method underestimates the low-count CI in the right panel.

![Confidence band comparison](https://github.com/usnistgov/FCpy/tree/main/doc/images/FCpy_execTime.png?raw=True)


## Benchmarks

Using ```%timeit``` magic command in IPython. All examples run on a Dell Latitude 5420 laptop with Intel i7-1185G7 processor @ 3.00GHz.

Zero counts, zero background (compare no correction vs. correction):

```python
%timeit FC.FC_poisson(0, 0, 100, conf= 0.95, useCorrection= False)
# [Out]: 1.96 ms ± 101 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
```

```python
%timeit FC.FC_poisson(0, 0, 100, conf= 0.95, useCorrection= True)
# [Out]: 3.57 ms ± 109 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

Zero counts, background rate = 0.1:

```python
%timeit FC.FC_poisson(0, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 3.74 ms ± 68.1 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

```python
%timeit FC.FC_poisson(0, 0.1, 100, conf= 0.95, useCorrection= True)
# [Out]: 9 ms ± 336 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

Ten counts, background rate = 0.1:

```python
%timeit FC.FC_poisson(10, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 7.22 ms ± 103 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

```python
%timeit FC.FC_poisson(10, 0.1, 100, conf= 0.95, useCorrection= True)
# [Out]: 103 ms ± 1.53 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

Above 20 counts, the Roe & Woodroofe (1999) correction becomes much more computationally expensive and likely uneccesary.

```python
%timeit FC.FC_poisson(50, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 29.3 ms ± 424 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

```python
%timeit FC.FC_poisson(50, 0.1, 100, conf= 0.95, useCorrection= True)
# [Out]: 2.14 s ± 78.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Even at a large number of counts where a Gaussian approximation would suffice (e.g., n0 = 1000), the calculation speed is not prohibitive:

```python
%timeit FC.FC_poisson(1000, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 1.76 s ± 3.83 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

The ```mumin``` and ```mumax``` kwargs can be used to limit the range of mu to explore (be careful!).

```python
%timeit FC.FC_poisson(1000, 0.1, 100, conf= 0.95, useCorrection= False, mumin=900, mumax=1100)
# [Out]: 1.19 s ± 8.64 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

The highly vectorized numpy implementation is very fast, and likely sufficient for most use cases. I made a numba implementation (not included in this package), however, the jit compile time far exceeded the time savings for one-off uses. Time savings were negligible when calculating up to 100 CIs. 

Benchmark comparison using the ```timeit``` module. Performance is favorable even compared to ROOT ```TFeldmanCousins()``` object called through ```PyROOT```.

![Benchmark comparison](https://github.com/usnistgov/FCpy/tree/main/doc/images/FCpy_execTime.png?raw=True)




## Installation

Requirements to install this Python module:

-   Python 3.8 or newer
-   numpy
-   scipy

(Optional)

-   dask

(Optional, for GUI):

-   PyQt5

(Optional, for Bayesian example):
-   pymc
-   arviz

Installing from source:

```shell
$ python setup.py build
$ pip install -e .
```

## Bayesian Example

The Bayesian example uses the PyMC library to simulate the posterior of the true counts (or count rate) of the signal of interest with a known background rate & uncertainty. The background rate prior is modeled with a TruncatedGaussian() prior with a lower limit of 0; the counts prior is Uniform() between 0 and 1M (the upper limit is arbitrary). The sum of the background and counts priors are used as the rate parameter of a Poission log-likelihood model, with observed n total counts. The CI is given as the highest density interval (HDI), the minimum width interval that encloses the given % of the posterior area.

Comparison of FC and Bayesian results:

For 0 counts in 1000 s with mean background rate 0.01 +/- 0.0005 (5% relative error):
Bayesian 95% HDI: 0.0003 3.097

The FC CI:
```python
FC.FC_poisson(0, 0.01, 1000, conf= 0.95, useCorrection= False)
# [Out]: array([0. , 1.2619])
```

```python
FC.FC_poisson(0, 0.01, 1000, conf= 0.95, useCorrection= True)
# [Out]: array([0. , 3.7642])
```

For 10 counts in 100 s with mean background rate 0.01 +/- 0.0005 (5% relative error):
Bayesian 95% HDI: 3.962 16.760

The FC CI:
```python
FC.FC_poisson(10, 0.01, 100, conf= 0.95, useCorrection= False)
# [Out]: array([3.752, 16.816])
```

```python
FC.FC_poisson(10, 0.01, 100, conf= 0.95, useCorrection= True)
# [Out]: array([3.752, 16.816])
```

## References

1) Brun R. and Rademakers F. (1997) ROOT — An object oriented data analysis framework. Nuclear Instruments and Methods in Physics Research A 389, 81-86.
2) Brun R., Rademakers F., Canal P., Naumann A., Couet O., Moneta L., Vassilev V., Linev S., Piparo D., Ganis G., Bellenot B., Guiraud E., Amadio G., Wverkerke, Mato P., TimurP, Tadel M., Wlav, Tejedor E., Blomer J., Gheata A., Hageboeck S., Roiser S., Marsupial, Wunsch S., Shadura O., Bose A., CristinaCristescu, Valls X. and Isemann R. (2019) root-project/root: v6.18/02. Zenodo 10.5281/zenodo.3895860.
3) Coakley K. J., Splett J. D. and Simons D. S. (2010) Frequentist coverage properties of uncertainty intervals for weak Poisson signals in the presence of background. Measurement Science and Technology 21.
4) Feldman G. J. and Cousins R. D. (1998) Unified approach to the classical statistical analysis of small signals. Physical Review D 57, 3873-3889.
5) Roe B. P. and Woodroofe M. B. (1999) Improved probability method for estimating signal in the presence of background. Physical Review D 60.
6) Roe B. P. and Woodroofe M. B. (2000) Setting confidence belts. Physical Review D 63.


## Contact Information
Evan Groopman <evan.groopman@nist.gov>

Surface and Trace Chemical Analysis Group

Materials Measurement Science Division (MMSD)

Material Measurement Laboratory (MML)

National Institute of Standards and Technology (NIST)
