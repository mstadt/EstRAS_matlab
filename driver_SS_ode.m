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
do_ACEi = false;
do_ARB = false;

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
                        'do_ACEi', do_ACEi,...
                        'do_ARB', do_ARB),...
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
SSdat = y(end,:)'; % steady state is end point

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



