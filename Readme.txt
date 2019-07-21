Curve fitting for 1D vector with any given model equations.

Run 'example.m' for an example
change your fitting model in the function 'mdl' in 'jcfitf.m' 

cite:
https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.6b05697
Jixin Chen, Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky, Christy F Landes, A Two-Step Method for smFRET Data Analysis, J. Phys. Chem. B, 2016, 120 (29), pp 7128–7132

A semi-exhaustive and brutal searching algorithm using the least-square method to find the best fit of a curve. The least-square mehtod can be changed easily to other residual treatment methods.

The fitting time of each iterature is scaled with the number of parameters and the accuracy of the fitting. Each parameter is searched one by one within its boundaries. Thus a P^N question is reduced to PN with a cost of losing space coverage.

The code can be changed to fit multiple curves globally by introduce different models or a model calls different parameters.
an example can be find at:

Juvinch R. Vicente, Ali Rafiei Miandashti, Kurt Sy Piecco, Joseph R. Pyle, Martin E. Kordesch, Jixin Chen*, Single-Particle Organolead Halide Perovskite Photoluminescence as a Probe for Surface Reaction Kinetics. ACS Applied Matierals & Interfaces, 2019, 11(19), 18034-18043.

by Jixin Chen @ Ohio University 2019
