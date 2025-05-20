function dydt = mod_eqns(t,y,params,varargin)
% Model equations for estrogen-RAS model

% estrogen
do_EST = 0; % change estrogen with time
do_EST_RAS = do_EST; % estrogen impact on RAS
EST_pct = 1; % percent of baseline estrogen

% RASi
do_ACEi = false;
do_ARB = false;

% infusion experiments
do_ANGII_inf = false;

for ii = 1:2:length(varargin)
    if strcmp(varargin{ii}, 'do_ANGII_inf')
        do_ANGII_inf = varargin{ii+1};
    elseif strcmp(varargin{ii}, 'do_EST')
        temp = varargin{ii+1};
        do_EST = temp(1);
        if do_EST == 1
            if length(temp) == 2
                do_EST_RAS = temp(2);
            else
                do_EST_RAS = do_EST;
            end
        elseif do_EST == 2
            if length(temp) == 3
                do_EST_RAS = temp(2);
                EST_pct = temp(3);
            else 
                fprintf('length of do_EST entry: %i \n', length(temp))
                error('do_EST must have entry for do_EST, do_EST_RAS, and EST_pct')
            end
        end
    elseif strcmp(varargin{ii}, 'do_ACEi')
        temp = varargin{ii+1};
        do_ACEi = temp(1);
        if do_ACEi
            age_ACEi = temp(2);
            pct_inhib_ACEi = temp(3);
        end
    elseif strcmp(varargin{ii}, 'do_ARB')
        temp = varargin{ii+1};
        do_ARB = temp(1);
        if do_ARB
            age_ARB = temp(2);
            pct_inhib_ARB = temp(3);
        end
    else
        fprintf('What is this varargin input? %s \n', varargin{ii})
        error('varagin input not available')
    end
end

%% Set variable names
% RAS
PRC       = y(1); % 
AGT       = y(2); % pmol/ml
AngI      = y(3); % fmol/ml = pmol/L
AngII     = y(4); % fmol/ml= pmol/L
Ang17     = y(5); % fmol/ml= pmol/L
AngIV     = y(6); % fmol/ml= pmol/L
AT1R      = y(7); % fmol/ml= pmol/L
AT2R      = y(8); % fmol/ml= pmol/L

%% Set parameter names
% NOTE: update updatepars.m to get parameter list from set_params()
h_renin = params(1);
h_AGT = params(2);
k_AGT = params(3);
X_PRCPRA = params(4);
KAGT0 = params(5);
Nrs = params(6);
A = params(7);
B = params(8);
AT1R0 = params(9);
c_ACE = params(10);
c_Chym = params(11);
c_NEP = params(12);
c_ACE2 = params(13);
c_IIIV = params(14);
c_AT1R = params(15);
c_AT2R = params(16);
h_AngI = params(17);
h_AngII = params(18);
h_Ang17 = params(19);
h_AngIV = params(20);
h_AT1R = params(21);
h_AT2R = params(22);
eRen = params(23);
eAGT = params(24);
eACE = params(25);
eAT1R = params(26);
eAT2R = params(27);


%% Model equations
dydt = zeros(length(y),1);

% Estrogen
if do_EST == 1
    EST = get_estrogen(t); 
elseif do_EST == 0
    EST = 1;
elseif do_EST == 2
    EST = EST_pct;
else
    fprintf('WARNING: do_EST not set up correctly \n')
end



%% RAS model
    
% Do estrogen impacts on RAS
if do_EST_RAS
    % Estrogen inhibits renin secretion (inhibition)
    etaEST_Ren_max = (eRen + 1)/eRen; 
    etaEST_Ren =  etaEST_Ren_max*g_minus(EST/eRen);
    % Estrogen impact on AGT (activation)
    etaEST_AGTmax = eAGT + 1;
    etaEST_AGT = etaEST_AGTmax * g_plus(EST/eAGT);
    % Estrogen impact on ACE (inhibition)
    etaEST_ACEmax = (eACE + 1)/eACE;
    etaEST_ACE = etaEST_ACEmax * g_minus(EST/eACE);
    % Estrogen impact on AT1R receptor binding (inhibition)
    etaEST_AT1Rmax = (eAT1R + 1)/eAT1R;
    etaEST_AT1R = etaEST_AT1Rmax * g_minus(EST/eAT1R);
    % Estrogen impact on AT2R receptor binding (activation)
    etaEST_AT2Rmax = (eAT2R + 1);
    etaEST_AT2R = etaEST_AT2Rmax * g_plus(EST/eAT2R);
    %etaEST_AT2R = 1;
else
    etaEST_AGT = 1;
    etaEST_Ren = 1;
    etaEST_AT1R = 1;
    etaEST_ACE = 1;
    etaEST_AT2R = 1;
end

eta_exp = A - B * log10(AT1R/AT1R0);
eta_AT1 = 10.^eta_exp;


R_sec = Nrs * eta_AT1 * etaEST_Ren;
% d(PRC)/dt
dydt(1) = R_sec - (log(2)/h_renin) * PRC;


PRA = PRC * X_PRCPRA * 2 * (AGT/(AGT + KAGT0));
% d(AGT)/dt
dydt(2) = etaEST_AGT * k_AGT - PRA - (log(2)/h_AGT)*AGT;

if do_ACEi
    t_yrs = t/(365*24);
    if t_yrs > age_ACEi
        if t_yrs < age_ACEi + 1
            % gradual decrease for numerics
            pct_ACEi = (t_yrs - age_ACEi)*pct_inhib_ACEi;
        else
            pct_ACEi = pct_inhib_ACEi;
        end
    else
        pct_ACEi = 0;
    end
else
    pct_ACEi = 0;
end
AngItoAngII = (etaEST_ACE*(1 - pct_ACEi)*c_ACE + c_Chym)*AngI;
AngItoAng17 = c_NEP * AngI;
% d(AngI)/dt
dydt(3) = PRA - AngItoAngII - AngItoAng17 - (log(2)/h_AngI)*AngI;

if do_ANGII_inf
    kinf0 = 130;
    % Graded infusion (Grant 1992)
    if t < 0
        k_inf2 = 0;
    elseif t < 1
        if t < (20/60)
            k_inf2 = kinf0;
        elseif t < (40/60)
            k_inf2 = kinf0*3;
        else
            k_inf2 = kinf0 * 10;
        end
    else
        k_inf2 = 0;
    end
else
    k_inf2 = 0;
end

if do_ARB
    t_yrs = t/(365*24);
    if t_yrs > age_ARB
        if t_yrs < age_ARB + 1
            % gradual decrease for numerics
            pct_ARB = (t_yrs - age_ARB)*pct_inhib_ARB; 
        else
            pct_ARB = pct_inhib_ARB;
        end
    else
        pct_ARB = 0;
    end
else
    pct_ARB = 0;
end
% d(AngII)/dt
AngIItoAT1R = etaEST_AT1R*(1 - pct_ARB)*c_AT1R * AngII;
dydt(4) = AngItoAngII ...
           - (c_ACE2 + c_IIIV +...
           + etaEST_AT2R*c_AT2R)*AngII ...
           - AngIItoAT1R...
           - (log(2)/h_AngII)*AngII ...
       + k_inf2;

% d(AngI7)/dt
dydt(5) = AngItoAng17 + c_ACE2 * AngII...
    - log(2)/h_Ang17 * Ang17;

% d(AngIV)/dt
dydt(6) = c_IIIV * AngII - log(2)/h_AngIV * AngIV;

% d(AT1R)/dt
dydt(7) = AngIItoAT1R - log(2)/h_AT1R * AT1R;

% d(AT2R)/dt
dydt(8) = etaEST_AT2R*c_AT2R * AngII - log(2)/h_AT2R * AT2R;

end % mod_eqns