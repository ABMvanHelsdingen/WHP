# WHP

## Cues
 
Contains files relating to the fitting of the whale cue data and the results. 
Note that the cue data itself is not public.


Bayesian.R: code to run Markov Chain Monte Carlo (MCMC) chain in the R package NIMBLE. Fits cue data to first model.<br />
Bayesian-Mixture.R: as above but for second, mixture model.<br />
Frequentist.R: code to fit cue data to first model by MLE in the R package TMB.<br />
whales_cues.cpp: TMB C++ template for fitting the first model, with baseline given by Equation 13, to the whale cue data by MLE.<br />

## Figures

Contains plots of the simulation study results

WIHP-<param>byK.jpeg: Histograms of the bias in <param> grouped by the true value of k for the first simulation study.<br />
MWIHP-<param>byK1.jpeg: Histograms of the bias in <param> grouped by the true value of k_1 for the second simulation study.<br />
MWIHP-K2byK2.jpeg: Histograms of the bias in K2 grouped by the true value of k_2 for the second simulation study.<br />


## SimResults

Contains the raw results of the two simulation studies.

WIHPsims_i_.csv (where i = 1:64) contain simulation results for the first simulation study.<br />
MWIHPsims_i_.csv (where i = 1:192) contain simulation results for the second simulation study, being the output of MWIHPsimulations.R<br /> 
They contain the  matrix output from WIHPsimulations.R and MWIHPsimulations.R<br />

## Simulations

Contains files for running the two simulation studies.

WHPsimulations.R: simulates the first model when the baseline rate mu(t) is constant (not used in manuscript).<br />
WIHPsimulations.R: simulates the first model when the baseline rate is inhomogeneous.<br />
MWIHPsimulations.R: simulates the mixture model.<br />

Simulate.R: contains functions called by the three above R scripts that implement the algorithms described in Appendix C. 

weibull_hawkes.cpp: TMB template for fitting the first model (with constant baseline) by MLE.<br />
weibull_hawkes_inhomog.cpp: TMB template for fitting the first model with a baseline function used in the simulation studies.<br /> 
mix_weibull_hawkes_inhomog.cpp: As above, but for the mixture model.<br />
