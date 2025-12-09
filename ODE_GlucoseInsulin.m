function dydt = ODE_GlucoseInsulin(t, y, params, funs)
% glucose and insulin on short time scale, including gluconeogenesis,
% assuming fixed number of insulin-producing beta-cells

% ----- dynamic variables
G = y(1); % glucose
I = y(2); % insulin


% ----- parameters
diet = params(1); % rate of glucose intake from diet
dG0 = params(2); % passive glucose disposal (e.g. brain, liver)
dGI = params(3); % maximal rate of insulin-stimulated glucose disposal
KI = params(4); % insulin resistance; insulin for half-maximal glucose disposal
pIB = params(5); % maximal rate of insulin secretion, assuming fixed number of beta cells
KG = params(6); % glucose level for half-maximal insulin secretion
dI = params(7); % insulin decay rate
gcn = params(8); % maximal rate of gluconeogenesis
KIg = params(9); % insulin for half-maximal gluconeogenesis
%KIg = KI + 10; % to make gluconeogenesis susceptible to insulin resistance when perturbing KI
nIg = params(10); % Hill coeff for gluconeogenesis
nI = params(11); % Hill coeff for insulin-stimulated glucose disposal
nG = params(12); % Hill coeff for glucose-stimulated insulin secretion 


diet_ = funs.diet(t,diet); % for oscillations simulating mealtimes


% -------- ODEs -------------------------------------------------------

% dG/dt: glucose
dydt(1) = diet_ + gcn*(KIg^nIg)/(KIg^nIg + I^nIg) - G*(dG0 + dGI*(I^nI)/(I^nI + KI^nI));
%        eating         gluconeogenesis        passive & insulin-stim glucose disposal

% dI/dt: insulin
dydt(2) = pIB*(G^nG)/(G^nG + KG^nG) - dI*I;
%                 GSIS             decay

dydt = dydt';
end


% FASTING blood glucose: < 5.6 mmol/L normally (and >4.0), < 7.0 mmol/L diabetes
% Plasma glucose peaks 1-2 hours after meal
% TWO HOURS AFTER meal, glucose < 7.8 mmol/L normally, 7.8-11.0 mmol/L prediabetes, and > 11.1 mmol/L diabetes
% POSTPRANDIAL INSULIN (microU/L) range 30-230 (@30 min), 18-276 (@1hr), 16-166 (@2hr)
