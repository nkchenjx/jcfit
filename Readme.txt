Curve fitting for 1D vector with any given model equations.

Developed and tested on software versions
MATLAB 2014b
Windows 10 

1. in folder jcfit_L1
    run jcfit_L1_test.m

2. in folder jcfit_L2
    run jcfittest.m

Both are the same code using different penalty terms:
1. L1 minimize the absolute sum of residual, i.e. least absolute sum of variation 
2. L2 minimize sum of square of residual, i.e. least square regression

They both have different strengths for different problems. 
Least square regression jcfit_L2 is >10 times slower than the Levenberg-Marquardt method for simple models such as exponential decays but is extremely easy to understand and to extend for complicated systems such as global fittings, fitting of multidimensional data, or systems with combined models.

Change your fitting model in the function 'mdl' in 'jcfitf.m' (moved to older versions 11/2021) if you use custom functions for a complicated model or global fitting.

cite:
https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.6b05697
Jixin Chen, Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky, Christy F Landes, A Two-Step Method for smFRET Data Analysis, J. Phys. Chem. B, 2016, 120 (29), pp 7128–7132

A semi-exhaustive and brutal searching algorithm using the least-square method to find the best fit of a curve. The least-square method can be changed easily to other residual treatment methods.

The fitting time of each iteration is scaled with the number of parameters of the fitting, the search speed is not very sensitive to the accuracy because of the exponential searching spacing design (not set as an option now and has to be changed manually inside the function). Each parameter is searched one by one within its boundaries. Thus a P^N question is reduced to PN with a cost of losing space coverage.

The code can be changed to fit multiple curves globally by introducing different models for different curves using the same parameters.

An example can be found at:
Juvinch R. Vicente, Ali Rafiei Miandashti, Kurt Sy Piecco, Joseph R. Pyle, Martin E. Kordesch, Jixin Chen*, Single-Particle Organolead Halide Perovskite Photoluminescence as a Probe for Surface Reaction Kinetics. ACS Applied Matierals & Interfaces, 2019, 11(19), 18034-18043.

by Jixin Chen @ Ohio University 2021
