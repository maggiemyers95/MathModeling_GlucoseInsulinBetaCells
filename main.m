clear all; close all;

%% ------ for beta-cell model single-param sensitivity analysis

% for each dynamic variable, define y-axis label, plot title, y-axis scale, y-axis limits (not used), and x-axis labels
plotAttributes = {
    struct('ylabel','Glucose (mmol/L)','title','Plasma Glucose Concentration','yscale','linear','ylim',[3,12],'xlabel','Time (days)'),
    struct('ylabel','Insulin (\muU/mL)', 'title','Plasma Insulin Concentration','yscale','linear','ylim',[3,200],'xlabel','Time (days)'),
    struct('ylabel','Log_{10} \beta_1-cells','title','Immature \beta-cell Population','yscale','log','ylim',[8,10],'xlabel','Time (days)'),
    struct('ylabel','Log_{10} \beta_2-cells','title','Mature \beta-cell Population','yscale','log','ylim',[8,10],'xlabel','Time (days)'),
    struct('ylabel','Log_{10} Total \beta-cells','title','Total \beta-cell Population','yscale','log','ylim',[8,10],'xlabel','Time (days)'),
    struct('ylabel','Percent \beta_2-cells','title','Proportion Mature \beta-cells','yscale','linear','ylim',[0,100],'xlabel','Time (days)')
    };

model_structure = 'pG2p'; % 'none'
funs = get_regulatory_functions(model_structure); % specifies beta-cell regulation

[params, y0, bounds, paramNames, ode_function] = get_parameters('FasterBetaTurnover'); nParams = length(paramNames);
funs.diet = @(t,p) p; % assume constant influx of glucose from diet, ideal for longer timescales
%funs.diet = @(t,p) 15.6755*p.*heaviside((sin(4*pi.*t) - 0.98));  % 45-minute meal twice a day
funs.KI = @(t,KI) KI; % constant insulin sensitivity
timespan = 0:0.1:365; % one year
timespan = 0:1:5*365; 
nVals = 5; % number of values to try for each perturbed parameter
which_params = [4,11,6]; % index of ordered list pars & paramNames
which_YtoPlot = [1,2;5,6]; % corresponds to each struct in plotAttributes


SingleParameterSensitivityAnalysis(params, bounds, nVals, which_params, which_YtoPlot, ode_function, funs, y0, timespan, paramNames, plotAttributes,'northeast')

which_YtoPlot = [1,3,4];

SingleParameterSensitivityAnalysis(params, bounds, nVals, which_params, which_YtoPlot, ode_function, funs, y0, timespan, paramNames, plotAttributes,'northeast')


%% ---- glucose and insulin only, on short timescale with oscillating diet to simulate meals


% FASTING blood glucose: < 5.6 mmol/L normally (and >4.0), < 7.0 mmol/L diabetes
% Plasma glucose peaks 1-2 hours after meal
% TWO HOURS AFTER meal, glucose < 7.8 mmol/L normally, 7.8-11.0 mmol/L prediabetes, and > 11.1 mmol/L diabetes
% POSTPRANDIAL INSULIN (microU/L) range 30-230 (@30 min), 18-276 (@1hr), 16-166 (@2hr)

close all;

[params, y0, bounds, paramNames, ode_function] = get_parameters('GluIns'); nParams = length(paramNames);
which_params = 1:nParams;
which_YtoPlot = [1,2];
nVals = 5;
timespan = 0:0.01:1; % one day
funs.diet = @(t,p) 15.6755*p.*heaviside((sin(4*pi.*t) - 0.98));  % 45-minute meal twice a day
funs.diet = @(t,p) p.*heaviside((sin(4*pi.*t) - 0.98));  % 45-minute meal twice a day

SingleParameterSensitivityAnalysis(params, bounds, nVals, which_params, which_YtoPlot, ode_function, funs, y0, timespan, paramNames, plotAttributes,'northeast')



%% ------- look at oscillating functions ------

%funs.diet = @(t,p) 15.6755*p.*heaviside((sin(4*pi.*t) - 0.98));  % 45-minute meal twice a day
figure; 
plot(timespan, funs.diet(timespan, params(1)),'k','Linewidth',2 ); 
set(gca, 'YLim', [-1, 1.1*max(funs.diet(timespan, params(1)))], 'color','w');
xlabel('Time (days)');
ylabel('Dietary Glucose Influx');
title('Simulating Mealtimes');

