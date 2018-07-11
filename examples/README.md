# Description of the example files

In this folder are locate a series of example files describing how to use the fitting library in sort-of real-life scenarios (some more than others). All the data to fit is generated inside the script themselves.

The examples and their content are listed here below, in order of complexity.

## Examples for `Fit`

- `ex_3expDecay.m`: fitting of a decay curve which is the sum of three exponential components. The data is supposes to come from some photon-counting experiment, so it is sampled from a Poisson distribution.  
The use of lower and upper bounds for the fit parameters is demonstrated.
- `ex_3expDecayIRF.m`: a variation of the above script, here the data is supposed to be measured on a system with a finite Impulse Response Function (IRF). The ability of `Fit` to work transparently while taking into account the IRF is demonstrated.

## Examples for `GlobalFitSimple`

- `ex_GFS_linear_fit.m`: global fit of a set of data following a linear dependence on the independent variable. The slope is equal for all data sets, while the vertical shift is different for each data set.  
The basic use of `GlobalFitSimple` is demonstrated together with generation of data and visual comparison of the data and fit.
- `ex_GFS_exp_decays.m`: we have a series of decay curves (sum of two exponentials) that are characterized by the same pair of decay times, while the amplitude of the two components change for different datasets.  
In this case the visual comparison is made both between the data and fit and on the parameters ("true" vs. obtained from the fit).
- `ex_GFS_Gaussians.m`: we fit a series of Gaussians whose position and width are fixed, while the amplitudes change among the datasets.
- `ex_GFS_Gaussians2.m`: this time the amplitude and width are kept constant, while the position varies with dataset.  
Here the use of lower/upper bound as well as th possibility to fix parameters is introduced.
- `ex_GFS_multipressure_Voigt.m`: the global fitting library is used to fit data simulated from spectral absorption lines with a Voigt profile. A linear pressure dependence of the position and width is supposed and fitted.
