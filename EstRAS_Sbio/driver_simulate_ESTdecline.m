% Simulate model with estrogen decline
clearvars; % clear

% Load simbiology model
model = copyobj(sbioloadproject("RAS_EST_v1.sbproj").m1);

sbioaccelerate(model) % accelerate model
getequations(model) % print model equations

% Change stop time
tf = 80*24*365; % 80 years
configset = getconfigset(model);
set(configset,'StopTime', tf);

% Get the solver options 
solverOptions = get(configset, 'SolverOptions');
% Set the maximum time step 
set(solverOptions, 'MaxStep', 0.5*365*24); % at least 2 time points per year

% turn on estrogen decline
rule = model.Rules(1);
rule.Active = true;

% Simulate
[time,x,names] = sbiosimulate(model);

% Time in years
t = time/(365*24);

%% Plot results
lw = 4;
fsize = 18;
xlab = 't (years)';
xminmax = [20,80];

figure(1);
clf;
tiledlayout(2,3);

% PRC
nexttile(1);
id = 1;
plot(t, x(:,id), 'linewidth', lw)
xlabel(xlab)
ylabel(names{id})
xlim(xminmax)
set(gca,'fontsize',fsize)
grid on

% AGT
nexttile(2);
id = 2;
plot(t, x(:,id)/1000, 'linewidth', lw)
xlabel(xlab)
xlim(xminmax)
ylabel(names{id})
set(gca,'fontsize',fsize)
grid on

% Ang I
nexttile(3);
id = 3;
plot(t, x(:,id), 'linewidth', lw)
xlabel(xlab)
xlim(xminmax)
ylabel(names{id})
set(gca,'fontsize',fsize)
grid on


% Ang II
nexttile(4);
id = 4;
plot(t, x(:,id), 'linewidth', lw)
xlabel(xlab)
xlim(xminmax)
ylabel(names{id})
set(gca,'fontsize',fsize)
grid on

% AT1R
nexttile(5);
id = 7;
plot(t, x(:,id), 'linewidth', lw)
xlabel(xlab)
ylabel(names{id})
xlim(xminmax)
set(gca,'fontsize',fsize)
grid on

% AT2R
nexttile(1);
id = 8;
plot(t, x(:,id), 'linewidth', lw)
xlabel(xlab)
ylabel(names{id})
xlim(xminmax)
set(gca,'fontsize',fsize)
grid on

