% Compute model steady state
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

% Set estrogen level

% turn off estrogen decline
rule = model.Rules(1);
rule.Active = false;

% set estrogen level at fixed level
EST_pct = 1;
param = sbioselect(model, "Type", "parameter","Name","EST");
param.Value = EST_pct; % set EST to fixed value


% compute steady state
[success, variant_out, mod_out, exitInfo] = sbiosteadystate(model);

disp(exitInfo)
disp(mod_out.Species)

% get species information
speciesList = sbioselect(mod_out, 'Type', 'Species');
speciesNames = {speciesList.Name};
SS_values = [speciesList.InitialAmount]; % Steady-state values

% Make a bar plot with the values for the species
fsize = 18;
figure(1);
clf;
tiledlayout(2,3);
for ii = [1,2,3,4,7,8]
    nexttile;
    bar(SS_values(ii))
    set(gca,'XTickLabel', speciesNames{ii},...
            'fontsize', fsize)
    ylabel('Steady-state concentration')
    grid on
end