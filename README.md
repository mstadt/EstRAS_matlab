# Estrogen-RAS model (MATLAB)
This code is for the May 2025 University of Waterloo Applied Mathematics Hackathon project. This code runs the estrogen-RAS model as presented in [Stadt and Layton, "A physiology-based mathematical model of the renin-angiotensin system, bone remodeling, and calcium homeostasis: Effects of estrogen and renin-angiotensin system inhibitors", submitted](https://www.biorxiv.org/content/10.1101/2025.05.01.651663v2.abstract)

## Script files
**driver_SS.m** Run this file to compute the steady state of the model using the function fsolve

**driver_ESTdecline.m** Run this file to conduct simulations of the RAS during age-related estrogen decline with and without RASi

## Functions
**mod_eqns** contains model equations

**set_params** sets up a struct with parameter names and values

**pars2vector** converts output of set_params() to a vector for use in MATLAB ODE solvers

**get_estrogen** age-related estrogen decline function

**updatepars** run this script to get a .txt file with the parameter names in the correct order from set_params(). This may be used to update mod_eqns when needing to change parameters
