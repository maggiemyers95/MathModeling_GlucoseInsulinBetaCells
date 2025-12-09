function [params, ICs, bounds, paramNames, ode_fun] = get_parameters(varargin)
% Returns parameter data to initialize ODE model 

% INPUT: Varargin (optional) can be used to retrieve variables corresponding
% to a different ODE model

% OUTPUT:
%     params: base-case parameters
%     ICs: initial conditions of each dynamic variable in ode_fun's equations
%     bounds: minimum and maximum parameter values
%     paramNames: cell of strings specifying each parameter name, using LaTeX math syntax
%     ode_fun: function handle referencing the file containing the system of ODEs


%--parameters format like [base-case, minimum, maximum] 
diet = [50, 30, 300]; % rate of glucose intake from diet
KI = [15, 15, 60]; % insulin resistance; insulin level for half-maximal glucose disposal
pI = [4e-6, 5e-7, 1e-5]; % maximal rate of insulin secretion per (mature) beta-cell
KG = [10, 7, 16]; % glucose for half-maximal insulin secretion 
dI = [20, 15, 30]; % insulin decay rate
dG0 = [2.5, 0.1, 7]; % rate of passive glucose disposal
dGI = [22.5, 10, 30]; % maximal rate of insulin-stimulated glucose disposal
gcn = [0, 0, 60]; % maximal rate of gluconeogenesis
KIg = [35, 15, 80]; % insulin for half-maximal gluconeogenesis

neo = [0, 0, 100]; % rate of beta-1 neogenesis
dB1 = [1/(30*365), 1/(40*365), 1/(3*365)]; % rate of beta-1 cell death
v = dB1(1).*[23, 5, 50]; % maximal rate of beta-1 mitosis... 23*dB1
p = [0.51, 0.35, 0.65]; % maximal probability of beta-1 self-renewal
q = [0.7, 0.4, 1]; % rate of beta-2 dedifferentiation (defined relative to dB2)
dB2 = [6, 5, 20]; % how much faster beta-2 cells die compared to beta-1 cells
m = [0, 0, 1]; % rate of mitosis-independent beta maturation, relative to mitosis rate
KGfb = [8, 5, 13]; % glucose for half-maximal regulation of beta-cells
KIfb = [15, 5, 25]; % insulin for half-maximal regulation of beta-cells
KB = [1e8, 3e7, 6e8]; % beta population for half-maximal regulation of beta-cells

%--parameter names (syntaxed for LaTeX interpreter)
paramNames = {'\alpha', 'd_{G_0}','d_{G_I}','K_I','p_I','K_G','d_I',...
    '\nu','v','p','q','d_{\beta_1}','d_{\beta_2}',...
    '\gamma','K_{I_\gamma}','m','K_{G^\ast}','K_{I^\ast}', 'K_{\beta^\ast}'};

% order of parsAll and paramNames must match ordered list of parameters defined in ode_fun 
parsAll = [diet;dG0;dGI;KI;pI;KG;dI;...
    neo;v;p;q;dB1;dB2;...
    gcn;KIg;m;KGfb;KIfb;KB]'; 


params = parsAll(1,:); % base-case parameters
bounds = parsAll(2:3,:); % min & max

ICs = [5,22,2.5e8, 7.5e8]; % initial conditions

ode_fun = @ODE_RegulateBetaCells; % function file containing model equations




% --------- 

if strcmp(varargin{1}, "GluIns") % model for only glucose and insulin (short timescales generally)
    
    %--parameters format like [base-case, minimum, maximum] 
    diet = [180, 50, 400]; % rate of glucose intake from diet
    KI = [15, 15, 60]; % insulin resistance; insulin level for half-maximal glucose disposal
    pIB = 5.*14*9e8.*[4e-6, 5e-7, 1e-5]; % maximal rate of insulin secretion per (mature) beta-cell
    KG = [5+10, 7, 16]; % glucose for half-maximal insulin secretion 
    dI = 7.5.*[20, 15, 30]; % insulin decay rate
    dG0 = [2.5, 0.1, 7]; % rate of passive glucose disposal
    dGI = [22.5, 10, 30]; % maximal rate of insulin-stimulated glucose disposal
    gcn = [170, 0, 250]; % maximal rate of gluconeogenesis
    KIg = [15, 10, 60]; % insulin for half-maximal gluconeogenesis
    nIg = [1,1,3]; % Hill coeff for gluconeogenesis, integer
    nI = [1,1,3]; % Hill coeff for insulin-stimulated glucose disposal, integer
    nG = [3,1,5]; % Hill coeff for glucose-stimulated insulin secretion

    % order of parsAll and paramNames must match ordered list of parameters defined in ode_fun 
    parsAll = [diet; dG0; dGI; KI; pIB; KG; dI; ...
        gcn; KIg; nIg; nI; nG]';

    % syntaxed for LaTeX interpreter
    paramNames = {'\alpha', 'd_{G_0}','d_{G_I}','K_I','p_I','K_G','d_I',...
    '\gamma','K_{I_\gamma}','n_\gamma','n_I','n_G'};

    % don't change these next 2 lines
    params = parsAll(1,:); % base-case parameters
    bounds = parsAll(2:3,:); % min & max

    ICs = [5,10]; % initial conditions
    
    ode_fun = @ODE_GlucoseInsulin; % function file containing model equations
end



if strcmp(varargin{1},"FasterBetaTurnover")
    
    %--parameters format like [base-case, minimum, maximum] 
    diet = [50, 30, 300]; % rate of glucose intake from diet
    KI = [15, 15, 60]; % insulin resistance; insulin level for half-maximal glucose disposal
    pI = [4e-6, 5e-7, 1e-5]; % maximal rate of insulin secretion per (mature) beta-cell
    KG = [10, 7, 16]; % glucose for half-maximal insulin secretion 
    dI = [20, 15, 30]; % insulin decay rate
    dG0 = [2.5, 0.1, 7]; % rate of passive glucose disposal
    dGI = [22.5, 10, 30]; % maximal rate of insulin-stimulated glucose disposal
    gcn = [0, 0, 60]; % maximal rate of gluconeogenesis
    KIg = [35, 15, 80]; % insulin for half-maximal gluconeogenesis
    
    neo = [0, 0, 100]; % rate of beta-1 neogenesis
    dB1 = [1/(2*365), 1/(10*365), 1/(1*365)]; % rate of beta-1 cell death
    v = dB1(1).*[23, 5, 50]; % maximal rate of beta-1 mitosis... 23*dB1
    p = [0.51, 0.35, 0.65]; % maximal probability of beta-1 self-renewal
    q = [0.7, 0.4, 1]; % rate of beta-2 dedifferentiation (defined relative to dB2)
    dB2 = [6, 5, 20]; % how much faster beta-2 cells die compared to beta-1 cells
    m = [0, 0, 1]; % rate of mitosis-independent beta maturation, relative to mitosis rate
    KGfb = [8, 5, 13]; % glucose for half-maximal regulation of beta-cells
    KIfb = [15, 5, 25]; % insulin for half-maximal regulation of beta-cells
    KB = [1e8, 3e7, 6e8]; % beta population for half-maximal regulation of beta-cells
    
    %--parameter names (syntaxed for LaTeX interpreter)
    paramNames = {'\alpha', 'd_{G_0}','d_{G_I}','K_I','p_I','K_G','d_I',...
        '\nu','v','p','q','d_{\beta_1}','d_{\beta_2}',...
        '\gamma','K_{I_\gamma}','m','K_{G^\ast}','K_{I^\ast}', 'K_{\beta^\ast}'};
    
    % order of parsAll and paramNames must match ordered list of parameters defined in ode_fun 
    parsAll = [diet;dG0;dGI;KI;pI;KG;dI;...
        neo;v;p;q;dB1;dB2;...
        gcn;KIg;m;KGfb;KIfb;KB]'; 
    
    
    params = parsAll(1,:); % base-case parameters
    bounds = parsAll(2:3,:); % min & max
    
    ICs = [5,22,2.5e8, 7.5e8]; % initial conditions
    
    ode_fun = @ODE_RegulateBetaCells; % function file containing model equations
end

end