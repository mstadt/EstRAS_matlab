% Run model simulation with age-related estrogen decline with 
% and without ACEi and ARB treatment
clearvars; % clear variables

%% load initial conditions
load("./IC/2025-05-20_ICfinal.mat","IC");

%% set parameters
p = set_params();
[params, ~] = pars2vector(p,0);

%% simulation settings
age_ACEi = 60; % years
pct_ACEi = 0.956;
age_ARB = 60; % years
pct_ARB = 0.9359;


%% Run simulation
% time span
t0 = 20*24*365; % Age 20 (hours)
tf = 80*24*365; % Age 80 (hours)
tspan = [t0,tf];

% ode solver options
opts_ode = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9, 'MaxStep', 10);

% No RASi simulation
fprintf("No RASi sim \n")
do_ACEi = 0;
do_ARB = 0;
do_EST = 1;
[t1, y1] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_EST', do_EST,...
                            'do_ACEi', [do_ACEi,age_ACEi,pct_ACEi],...
                            'do_ARB', [do_ARB,age_ARB,pct_ARB]),...
                            tspan, IC, opts_ode);

% ACEi simulation
fprintf("ACEi sim \n")
do_ACEi = 1;
do_ARB = 0;
do_EST = 1;
[t2, y2] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_EST', do_EST,...
                            'do_ACEi', [do_ACEi,age_ACEi,pct_ACEi],...
                            'do_ARB', [do_ARB,age_ARB,pct_ARB]),...
                            tspan, IC, opts_ode);

% ARB simulation
fprintf("ARB sim \n")
do_ACEi = 0;
do_ARB = 1;
do_EST = 1;
[t3, y3] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_EST', do_EST,...
                            'do_ACEi', [do_ACEi,age_ACEi,pct_ACEi],...
                            'do_ARB', [do_ARB,age_ARB,pct_ARB]),...
                            tspan, IC, opts_ode);

% convert time to years
t1 = t1./(365*24);
t2 = t2./(365*24);
t3 = t3./(365*24);

%% Plot results
% fig specs
cmap = parula(6);
c1 = cmap(1,:); ls1 = '-';
c2 = cmap(3,:); ls2 = '-.';
c3 = cmap(5,:); ls3 = '--';
lw = 5;
fsize = 18;
xminmax =[20,80];
labs = {'No RASi', 'ACEi', 'ARB'};
xlab = "t (years)";


fprintf("plotting results \n")
f = figure(1);
clf;
width = 1600;
height = 800;
f.Position = [100, 100, width, height];
tiledlayout(2,3);

% AGT
nexttile(1);
id = 2;
hold on;
plot(t1,y1(:,id)/1000, 'linewidth', lw, ...
                        'linestyle', ls1,...
                        'color', c1)
plot(t2,y2(:,id)/1000, 'linewidth', lw,...
                            'linestyle', ls2,...
                            'color',c2)
plot(t3,y3(:,id)/1000, 'linewidth',lw,...   
                            'linestyle', ls3,...
                            'color', c3)
grid on
ylabel('[AGT] (fmol/L)')
xlim(xminmax)
xlabel(xlab)
legend(labs)
set(gca,'fontsize',fsize)

% Renin (PRC)
nexttile(2);
id = 1;
hold on;
plot(t1,y1(:,id), 'linewidth', lw, ...
                        'linestyle', ls1,...
                        'color', c1)
plot(t2,y2(:,id), 'linewidth', lw,...
                            'linestyle', ls2,...
                            'color',c2)
plot(t3,y3(:,id), 'linewidth',lw,...   
                            'linestyle', ls3,...
                            'color', c3)
grid on
ylabel('[Renin] (mU/L)')
xlim(xminmax)
xlabel(xlab)
legend(labs, 'location','best')
set(gca,'fontsize',fsize)

% Ang I
nexttile(3);
id = 3;
hold on;
plot(t1,y1(:,id), 'linewidth', lw, ...
                        'linestyle', ls1,...
                        'color', c1)
plot(t2,y2(:,id), 'linewidth', lw,...
                            'linestyle', ls2,...
                            'color',c2)
plot(t3,y3(:,id), 'linewidth',lw,...   
                            'linestyle', ls3,...
                            'color', c3)
grid on
ylabel('[Ang I] (pmol/L)')
xlim(xminmax)
xlabel(xlab)
legend(labs, 'location','best')
set(gca,'fontsize',fsize)

% Ang II
nexttile(4);
id = 4;
hold on;
plot(t1,y1(:,id), 'linewidth', lw, ...
                        'linestyle', ls1,...
                        'color', c1)
plot(t2,y2(:,id), 'linewidth', lw,...
                            'linestyle', ls2,...
                            'color',c2)
plot(t3,y3(:,id), 'linewidth',lw,...   
                            'linestyle', ls3,...
                            'color', c3)
grid on
ylabel('[Ang II] (pmol/L)')
xlim(xminmax)
xlabel(xlab)
legend(labs, 'location','best')
set(gca,'fontsize',fsize)

% AT1R
nexttile(5);
id = 7;
hold on;
plot(t1,y1(:,id), 'linewidth', lw, ...
                        'linestyle', ls1,...
                        'color', c1)
plot(t2,y2(:,id), 'linewidth', lw,...
                            'linestyle', ls2,...
                            'color',c2)
plot(t3,y3(:,id), 'linewidth',lw,...   
                            'linestyle', ls3,...
                            'color', c3)
grid on
ylabel('[AT1R-bound Ang II] (pmol/L)')
xlim(xminmax)
xlabel(xlab)
legend(labs, 'location','best')
set(gca,'fontsize',fsize)

% AT2R
nexttile(6);
id = 8;
hold on;
plot(t1,y1(:,id), 'linewidth', lw, ...
                        'linestyle', ls1,...
                        'color', c1)
plot(t2,y2(:,id), 'linewidth', lw,...
                            'linestyle', ls2,...
                            'color',c2)
plot(t3,y3(:,id), 'linewidth',lw,...   
                            'linestyle', ls3,...
                            'color', c3)
grid on
ylabel('[AT2R-bound Ang II] (pmol/L)')
xlim(xminmax)
xlabel(xlab)
legend(labs, 'location','best')
set(gca,'fontsize',fsize)

AddLetters2Plots(figure(1), {'(A)', '(B)', '(C)', '(D)', '(E)','(F)'},...
                'HShift', -0.07, 'VShift', -0.06,...
                'fontsize', 20)


