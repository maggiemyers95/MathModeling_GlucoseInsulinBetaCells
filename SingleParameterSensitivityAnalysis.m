function SingleParameterSensitivityAnalysis(params, bounds, nVals, which_params, which_YtoPlot, ode_fun, funs, y0, timespan, paramNames, plotAttributes, varargin)
%-- Generate plots of ODE state variables over time, simulated for each color-coded parameter perturbation 
%-- INPUTS:  
%     params: (1 x nParams) specifying base-case parameter values
%     bounds: (2 x nParams) specifying [min; max] values of each parameter... empty [] is okay
%     nVals: integer number of values taken by each perturbed parameter
%     which_params: array of integers, each integer specifies the index of which parameter(s) to perturb
%     which_YtoPlot: array of integers specifying index of which model variable(s) to plot
%     ode_fun: function handle that specifies which ODE model to use
%     funs: cell containing optional input for ode_fun (allows for time-dependent parameter values)
%     y0: array of initial conditions
%     timespan: array of time points over which to simulate & plot model
%     paramNames: {1 x nParams} cell containing strings for each parameter name
%     plotAttributes: cell of structs that specify plot attributes (e.g. y-axis labels, titles) for each model variable
%     varargin: {1} legend position e.g. 'best' 'bestoutside'


%-- DETAILED SUMMARY: 
% for each parameter specified (by an integer) in which_params, perform nVal
% simulations over which that parameter's value is varied uniformly across given
% bounds. In plots, the black line corresponds to simulation with base-case
% parameters, and colorful lines correspond to simulations where a single
% parameter is perturbed, with the color choice corresponding to the
% parameter's value. 


%------------------------------------------------------------------------------------------------------------------------------------------------------
fontsize = 16; fontstyle = 'Helvetica';
nParams = length(paramNames); 
if timespan(end)<= 5 % shorter timescales typically use oscillating diet parameter to simulate mealtimes
    ode_options = odeset('MaxStep',0.015); % makes sure solver doesn't skip over mealtime
else 
    ode_options = [];
end

%-- base-case solution
sol = ode45(@(t,y)ode_fun(t,y,params,funs), timespan, y0, ode_options);
Y_base = deval(sol,timespan);
Y_base = Y_transform(Y_base);


%-- use parameter bounds if provided, otherwise vary from 50%-150% of parameter value
if isempty(bounds)
    for i = 1:nParams
        lowerBound(i) = 0.5*params(i); upperBound(i) = 1.5*params(i); 
    end
else
    lowerBound = bounds(1,:); upperBound = bounds(2,:);
end


%------------------------------------------------------------------------------------------------------------------------------------------------------
% --------- main loop for single-parameter sensitivity analysis ----------------------------------------------

for i = which_params % iterate each parameter to be perturbed
    perturbationValues = linspace(lowerBound(i),upperBound(i),nVals); % array of values taken by parameter i 
    infoString=sprintf('Perturbing %s = %g from [%g,  %g]',paramNames{i},params(i),lowerBound(i),upperBound(i));
    simParams = zeros(nVals,nParams); % initialize with zeros
    for iVal = 1:nVals % initialize with base parameters for each trial
        simParams(iVal,1:nParams) = params; 
    end
    simParams(:,i) = perturbationValues(1,:); % assign perturbation for each trial
    
    Yii = zeros(size(Y_base,1),length(timespan),nVals); % initialize solution arrays
    fprintf('Evaluating for %s...\t',paramNames{i});

    %-- iterate over each value for the parameter perturbed, simulate ODE, and record data for plotting
    for iVal = 1:nVals
        pars_i = simParams(iVal,:); % parameters for this simulation
        sol_i = ode45(@(t,y)ode_fun(t,y,pars_i,funs), timespan, y0, ode_options);
        Ys_i = deval(sol_i, timespan); % get solution y(t)
        Ys_i = Y_transform(Ys_i); % gets Btotal and percentB2, function allows flexibility to handle different ODE models
        Yii(:,:,iVal) = Ys_i; 
        fprintf('%g  ',iVal);
    end
    fprintf('\n'); 

    colorVec = 1:nVals; % each color corresponds to different parameter value
    %-- get y-axis bounds for each plot
    for k=1:size(Yii,1)
        Y_bounds(1,k) = min(min(Yii(k,:,:)));
        Y_bounds(2,k) = max(max(Yii(k,:,:)));
    end
    

    %---------- plot! ---------------------------------------------------------------------------------------------------------------------------------
    sp_dims = size(which_YtoPlot); % subplot dimensions (rows, columns)
    fig_width = 0.5*sp_dims(2); fig_height=0.6*sp_dims(1); % change multipliers if you want to
    figure('Name',sprintf('Sensitivity Analysis: %s',infoString)); 
    pos = get(gcf,'Position'); set(gcf,'Position', pos.*[1/(3*fig_width), 1/(3*fig_height), fig_width,fig_height]); hold on;
    counter=1; % for each subplot i.e. each element in which_toPlot
    for sp_i=1:sp_dims(1) % rows
        for sp_j=1:sp_dims(2) % columns
            k = which_YtoPlot(sp_i,sp_j); % gets index of ODE variable to plot i.e. 1=Glucose, 2=Insulin, 3=Beta1cells, ...
            subplot(sp_dims(1), sp_dims(2), counter); counter=counter+1;
            plot(timespan,Y_base(k,:), 'k:', 'LineWidth',2.5, 'DisplayName',sprintf('%g (base)',params(i))); hold on; % plot base case
            for j=1:nVals % for each parameter value
                plot(timespan,Yii(k,:,j), 'color', v2rgb( colorVec(j), colorVec),'linewidth',2, 'DisplayName',sprintf('%g',perturbationValues(1,j))); hold on;
            end
            title(plotAttributes{k}.title); xlabel(plotAttributes{k}.xlabel); ylabel(plotAttributes{k}.ylabel); 
            set(gca,'XLim',[0, timespan(end)], 'YScale', plotAttributes{k}.yscale); 
            set(gca, 'YLim',[0.9*Y_bounds(1,k), 1.1*Y_bounds(2,k)]); % plotAttributes{k}.ylim, 
            set(gca,'fontsize',fontsize,'fontname',fontstyle, 'color','w');
            hold off; 
        end
    end
    % legend attached to last subplot
    lgd = legend('FontSize',fontsize*0.8,'Location','best','color','w', 'textcolor','k', 'Interpreter','latex');
    title(lgd, sprintf('$%s$ value', paramNames{i}), 'FontSize',fontsize);
    if nargin >= 12
        lgd.Location = varargin{1};
    end

    drawnow
end


end