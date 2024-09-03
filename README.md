# Voids_RSD


This repository contains codes for the calculation of the real-space void properties such as density profiles (void-tracer correlation function), velocity profiles, and velocity dispersion profiles, and redshift space void properties such as the various multipole moments of the anisotropic void-tracer correlation function from the GLAM simuations as presented in 

C.Wilson and R. Bean (2022), https://arxiv.org/abs/2212.02569

From the above real space quantities, this code also estimates the redshift space density profile and its multipole moments through the use of the gaussian streaming model. Examples of how to use this code are provided in the exampleCalculations.py file in the analysis directory, which uses methods in the voids_MG_GSM_methods.py file also in the analysis directory. To reproduce plots from the main section of our paper, one may run voids_MG_GSM_plots.py

All code is run using highly standard python libraries. In addition to modules such as os, time, copy, etc... the following modules need to be available for use


numpy
matplotlib
scipy
seaborn

for any additional questions, please reach out to Christopher Wilson at cww78@cornell.edu



