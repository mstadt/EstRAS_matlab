function pars = set_params()
% set parameter values in a struct

%% RAS parameters
pars.h_renin = 12/60; % hours %12; % mins
pars.h_AGT = 10; % hours (Lo et al 2011)
pars.k_AGT = 610.39 * 60; % fmol/ml/hour % from female human get_pars in BP-regulation

pars.X_PRCPRA = 61/60 * 60; % fmol/hour/pg
pars.KAGT0 = 513525.02;%
pars.Nrs = 0.8*60; %
pars.A = 0.0102;  % A_AT1-renin, Hallow et al 2014
pars.B = 0.95;  % B_AT1-renin in BP-regulation code, Hallow et al 2014
pars.AT1R0 = 3.75; % female

pars.c_ACE = 1.4079 * 60; % (female BP-regulation) 1/hr 0.88492; % 1/min
pars.c_Chym = 0.1482 * 60; % (female BP-regulation) 1/hr % 0.09315; % 1/min
pars.c_NEP = 0.060759 * 60; %(female BP-regulation) 1/hr %0.038189; % 1/min

pars.c_ACE2 =  0.0037603 * 60; % (female, BP-regulation) 1/hr 0.0078009; % 1/min
pars.c_IIIV =  0.038644 * 60; % (female, BP-regulation) 1/hr 0.25056; % 1/min
pars.c_AT1R = 0.027089 * 60; % (female, BP-regulation) 1/hr 0.17008; % 1/min
pars.c_AT2R = 0.038699 * 60; % (female, BP-regulation) 1/hr %0.065667; % 1/min

pars.h_AngI = 0.5/60; % hours %0.5; % min
pars.h_AngII = 0.66/60; % hours %0.66; % min

pars.h_Ang17 = 30/60; % hours %30; % min
pars.h_AngIV = 0.5/60; % hours %0.5; % min
pars.h_AT1R  = 12/60; % hours 12; % min
pars.h_AT2R  = 12/60; % hours %12; % min



%% Estrogen impact on RAS
pars.eRen = 0.15; % estrogen impact on renin secretion
pars.eAGT = 0.15; % estrogen impact on AGT
pars.eACE = 5; % estrogen impact on ACE
pars.eAT1R = 0.5; % estrogen impact on AT1R
pars.eAT2R = 0.05;  % estrgen impact on AT2R
end