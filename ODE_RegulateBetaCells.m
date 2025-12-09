function dydt = ODE_RegulateBetaCells(t, y, params, funs)
% ODE model of glucose, insulin, beta-1 cells, and beta-2 cells


% ----- dynamic variables
G = y(1); % plasma glucose concentration
I = y(2); % plasma insulin concentration 
B1 = y(3); % population of immature beta-1 cells
B2 = y(4); % population of mature beta-2 cells


% ----- parameters
diet = params(1); % rate of glucose intake from diet
dG0 = params(2); % rate of passive glucose disposal
dGI = params(3); % maximum rate of insulin-stimulated glucose disposal
KI = params(4); % insulin resistance; insulin for half-maximal glucose disposal
pI2 = params(5);  %  maximum rate of insulin secretion by mature beta-2 cells
KG2_ = params(6);   % glucose for half-maximal insulin secretion by mature beta-2 cells
pI1 = 0.2*pI2; KG1 = 1.8*KG2_; % GSIS params of immature beta-1 cells relative to mature beta-2 cells
dI = params(7); % insulin decay rate
dB1 = params(12);  % rate of beta-1 cell death
dB2 = dB1*params(13); % rate of beta-2 cell death (defined relative to beta-1 death rate)
neo = params(8); % maximum rate of beta-1 cell neogenesis
v = params(9); % maximum rate of beta-1 cell division
p = params(10); % maximum probability of beta-1 cell self-renewal 
q = dB2*params(11); % maximum rate of beta-2 dedifferentiation
gcn = params(14); % maximum rate of gluconeogenesis
KIg = params(15); % insulin for half-maximal gluconeogenesis
%KIg = (KIg+KI)/2; % to make gluconeogenesis susceptible to changes in insulin resistance
m = v*params(16); % rate of mitosis-independent beta-1 cell maturation, relative to mitosis rate
KGfb = params(17); % glucose for half-maximal regulation of beta-cells
KIfb = params(18);  % insulin for half-maximal regulation of beta-cells
KB = params(19); % beta population for half-maximal regulation of beta-cells
n = 2; % Hill coeff for GSIS



% ----------- DIET & INSULIN RESISTANCE CHANGING OVER TIME? -------------
KI_ = funs.KI(t,KI); % insulin resistance
diet_ = funs.diet(t,diet); % diet
%diet__ = diet_funs{1}(t,diet_); % heaviside oscillations diet



% ------------ BETA-CELL REGULATION FUNCTIONS --------------------
neo_ = funs.neo(G,I,B1,B2,neo,KGfb,KIfb,KB); 
v_ = funs.v(G,I,B1,B2,v,KGfb,KIfb,KB);
p_ = funs.p(G,I,B1,B2,p,KGfb,KIfb,KB);
q_ = funs.q(G,I,B1,B2,q,KGfb,KIfb,KB);
m_ = funs.m(G,I,B1,B2,m,KGfb,KIfb,KB);



% ----- System of ODEs (ordinary differential equations) -------------

% dG/dt: glucose
dydt(1) = diet_  + gcn*(KIg)/(I + KIg) - G*(dG0 + dGI*I/(I + KI_));
%         eating   gluconeogenesis      passive & ins-stim glu disposal

% dI/dt: insulin
dydt(2) = pI1*B1*(G^n)/(G^n + KG1^n) + pI2*B2*(G^n)/(G^n + KG2_^n) - dI*I;
%         secretion by beta-1 cells    secretion by beta-2 cells     decay

% dB1/dt: Beta-1-cells (immature)
dydt(3) = neo_ + (2*p_ -1)*B1*v_   - m_*B1     - dB1*B1 + q_*B2;
%     neogenesis   cell division   maturation   death    dedifferentiation

% dB2/dt: Beta-2-cells (mature)
dydt(4) = 2*(1 - p_)*v_*B1    + m_*B1   - dB2*B2 - q_*B2;
%         maturation        maturation   death   dedifferentiation

dydt = dydt';

end