# FCpy

Calculate Feldman-Cousins confidence intervals with Python & numpy. 

Highly vectorized and fast implementation of Poisson and Gaussian confidence interval (CI) calculation algorithms described in Feldman & Cousins (1998) Phys. Rev. D, 57, 3873 https://doi.org/10.1103/PhysRevD.57.3873. Generally for Poisson processes in the presence of a known mean background count rate or Gaussian processes with a physical lower cutoff of 0. Note, these do not incorporate uncertainty on the mean background rate (see Bayesian example instead).

Includes optional low-count correction to Poisson CIs described by Roe & Woodroofe (1999) Phys. Rev. D, 60, 053009 https://doi.org/10.1103/PhysRevD.60.053009 (this is more computationally expensive, see Benchmarks). The original Feldman & Cousins formulation would result in tighter confidence intervals with increasing background levels, which does not make sense. This is especially problematic for the case of 0 total counts, where the number of background counts were clearly 0. If, for example, the average background rate were relatively large, then this measurement would not be representative of the true mean background - it would be a statistically unlikley, yet possible event. Therefore, the confidence on the non-background counts should not be reduced - this measurement was clearly an outlier.

Alternatively, a Bayesian example shows the use of PyMC to perform the analogous MCMC simulation. The benefits of the Bayesian solution are that it can incorporate uncertainties on the background rate and other instrumental sensitivity factors. The Bayesian solution also yields the CI of the parameter of interest (e.g., the percent probability that the number of counts lies within the CI), which is generally what the user actually wants. The frequentist interpretation is instead that upon repeating the same experiment X times, the result will fall within that range CI % of the repetitions.

Further reading: see Coakley et al. (2010) Meas. Sci. Tech. 21, 035102 https://doi.org/10.1088/0957-0233/21/3/035102 for a discussion of coverage properties between FC, Bayesian, and propagation of errors (POE) methods.


## Usage:

```shell
from FCpy.FCpy import FC
FC.FC_poisson(5, 0.01, 1000, conf= 0.95, useCorrection= False)
# [Out]: array([0. , 2.93706109])
```

```shell
FC.FC_poisson(5, 0.01, 1000, conf= 0.95, useCorrection= True)
# [Out]: array([0. , 4.76047702])
```

A convenience function calculates 68.3%, 95%, and 99.7% CIs:

```shell
FC.FC_poisson123(5, 0.01, 1000, useCorrection= False)
# [Out]: [(0.0, 0.27183253885023984), (0.0, 2.937061090328606), (0.0, 6.40224918079424)]
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


## Benchmarks

Using ```%timeit``` magic command in IPython. All examples run on a Dell Latitute 5420 laptop with Intel i7-1185G7 processor @ 3.00GHz.

Zero counts, zero background:

```shell
%timeit FC.FC_poisson(0, 0, 100, conf= 0.95, useCorrection= False)
# [Out]: 1.96 ms ± 101 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
```

```shell
%timeit FC.FC_poisson(0, 0, 100, conf= 0.95, useCorrection= False)
# [Out]: 3.57 ms ± 109 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

Zero counts, background rate = 0.1:

```shell
%timeit FC.FC_poisson(0, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 3.74 ms ± 68.1 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

```shell
%timeit FC.FC_poisson(0, 0.1, 100, conf= 0.95, useCorrection= True)
# [Out]: 9 ms ± 336 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

Ten counts, background rate = 0.1:

```shell
%timeit FC.FC_poisson(10, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 7.22 ms ± 103 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

```shell
%timeit FC.FC_poisson(10, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 103 ms ± 1.53 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

Above 20 counts, the Roe & Woodroofe (1999) correction becomes much more computationally expensive and likely uneccesary.

```shell
%timeit FC.FC_poisson(50, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 29.3 ms ± 424 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

```shell
%timeit FC.FC_poisson(50, 0.1, 100, conf= 0.95, useCorrection= True)
# [Out]: 2.14 s ± 78.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Even at a large number of counts where a Gaussian approximation would suffice (e.g., 1000), the calculation speed is not prohibitive:

```shell
%timeit FC.FC_poisson(1000, 0.1, 100, conf= 0.95, useCorrection= False)
# [Out]: 1.48 s ± 72.8 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```


## Installation

Requirements to install this Python module:

-   Python 3.8 or newer
-   numpy
-   scipy

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

```shell
FC.FC_poisson(0, 0.01, 1000, conf= 0.95, useCorrection= False)
# [Out]: array([0. , 1.2619])
```

```shell
FC.FC_poisson(0, 0.01, 1000, conf= 0.95, useCorrection= True)
# [Out]: array([0. , 3.7642])
```

For 10 counts in 100 s with mean background rate 0.01 +/- 0.0005 (5% relative error):
Bayesian 95% HDI: 3.962 16.76

```shell
FC.FC_poisson(10, 0.01, 100, conf= 0.95, useCorrection= False)
# [Out]: array([3.752, 16.816])
```

```shell
FC.FC_poisson(10, 0.01, 100, conf= 0.95, useCorrection= True)
# [Out]: array([3.752, 16.816])
```


## Contact Information
Evan Groopman <evan.groopman@nist.gov>

Surface and Trace Chemical Analysis Group

Materials Measurement Science Division (MMSD)

Material Measurement Laboratory (MML)

National Institute of Standards and Technology (NIST)
