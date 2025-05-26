% Get steady state of model for given parameter set
clearvars; % clear variables

% set parameters
p = set_params();

% change params here

[params, parnames] = pars2vector(p, 0);

%% Optional model inputs
do_EST = 2; % do estrogen effects - 0 - no estrogen, 1: age-related estrogen decline, 2: fix estrogen at EST_pct
do_EST_RAS = 1; % do estrogen effects on RAS
EST_pct = 1; % percent of baseline estrogen
do_ACEi = 1;
pct_inhib_ACEi = 0.9;
age_ACEi = 0;
do_ARB = 1;
pct_inhib_ARB = 0.9;
age_ARB = 0;

%% Run simulation from initial condition
% set initial conditions
load("./IC/2025-05-20_ICfinal.mat", "IC");

fprintf("solving ODEs \n")

% timespan
t0 = 0;
tf = 1e5;

tspan = [t0, tf];
opts_ode = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9, 'MaxStep', 10);

[t,y] = ode15s(@(t,y) mod_eqns(t,y,params,...
                        'do_EST', [do_EST, do_EST_RAS, EST_pct],...
                        'do_ACEi', [do_ACEi, age_ACEi, pct_inhib_ACEi],...
                        'do_ARB', [do_ARB, age_ARB, pct_inhib_ARB]),...
                        tspan, IC, opts_ode);

% Plot figure
% fig specs
lw = 4;
cmap = summer(5);
c1 = cmap(3,:);
fsize = 16;
xlab = 't';

figure(1);
clf;

nexttile;
id = 1;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('PRC')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 2;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('AGT')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 3;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('AngI')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 4;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('AngII')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 5;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('Ang17')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 6;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('AngIV')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 7;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('AT1R')
grid on
set(gca,'fontsize',fsize)

nexttile;
id = 8;
plot(t,y(:,id),'linewidth',lw)
ymin = min(y(:,id) - 1);
ymax = max(y(:,id) + 1);
ylim([ymin, ymax]);
xlabel(xlab)
ylabel('AT2R')
grid on
set(gca,'fontsize',fsize)

%% Find SS solution
% Use final point in ODE solution as initial guess
IG = y(end,:)';  % initial guess
opts_fsolve = optimoptions('fsolve','Display', 'none',...
                                'MaxFunEvals', 1e6,...
                                'MaxIter', 1e6,...
                                'FunctionTolerance', 1e-16);
[SSdat, residual, ...
    exitflag, output] = fsolve(@(y) mod_eqns(0, y, params,...
                                        "do_EST", [do_EST, do_EST_RAS, EST_pct],...
                                        "do_ACEi", [do_ACEi, age_ACEi, pct_inhib_ACEi],...
                                        "do_ARB", [do_ARB, age_ARB, pct_inhb_ARB]),...
                                        IG, opts_fsolve);

fprintf('maximum residual size: %0.1d \n', max(abs(residual)))
if exitflag ~= 1
    fprintf('WARNING: SS did not converge exitflag: %i \n', exitflag)
end

% Print out solution
if exitflag == 1
    fprintf('SS solution \n')
    fprintf('PRC:    %0.2f \n', SSdat(1))
    fprintf('AGT:    %0.2f \n', SSdat(2))
    fprintf('AngI:   %0.2f \n', SSdat(3))
    fprintf('AngII:  %0.2f \n', SSdat(4))
    fprintf('Ang17:  %0.2f \n', SSdat(5))
    fprintf('AngIV:  %0.2f \n', SSdat(6))
    fprintf('AT1R:   %0.2f \n', SSdat(7))
    fprintf('AT2R:   %0.2f \n', SSdat(8))
end

%% Bar plot with steady state results
f = figure(2);
clf;
width = 1600;
height = 800;
f.Position = [100, 100, width, height];
tiledlayout(2,3);

% AGT
nexttile(1);
id = 2;
hold on;
bar(SSdat(id)/1000)
xticks([])
grid on
ylabel('[AGT] (fmol/L)')
set(gca,'fontsize',fsize)

% Renin (PRC)
nexttile(2);
id = 1;
hold on;
bar(SSdat(id))
grid on
xticks([])
set(gca,'fontsize',fsize)
ylabel('[Renin] (mU/L)')


% Ang I
nexttile(3);
id = 3;
hold on;
bar(SSdat(id))
grid on
xticks([])
set(gca,'fontsize',fsize)
ylabel('[Ang I] (pmol/L)')


% Ang II
nexttile(4);
id = 4;
hold on;
bar(SSdat(id))
grid on
xticks([])
set(gca,'fontsize',fsize)
ylabel('[Ang II] (pmol/L)')


% AT1R
nexttile(5);
id = 7;
hold on;
bar(SSdat(id))
grid on
xticks([])
set(gca,'fontsize',fsize)
ylabel('[AT1R-bound Ang II] (pmol/L)')

% AT2R
nexttile(6);
id = 8;
hold on;
bar(SSdat(id))
grid on
xticks([])
set(gca,'fontsize',fsize)
ylabel('[AT2R-bound Ang II] (pmol/L)')


AddLetters2Plots(figure(2), {'(A)', '(B)', '(C)', '(D)', '(E)','(F)'},...
                'HShift', -0.07, 'VShift', -0.07,...
                'fontsize', 20)



