function funs = get_regulatory_functions(model_structure)
% Use a Hill (or Michaelis-Menten) function to model regulatory control of
% a process/parameter. The Hill function takes a value in [0,1] and 
% represents the percent activation of the process being regulated.

% INPUTS:
%      model_structure: string specifying the type of beta-cell regulation... 
%          * 1st character: the parameter to regulate, choose from {p,v,q,m,n}
%          * 2nd character: the dynamic model variable controlling/applying the regulation, choose from {G,I,B1,B2}
%          * second-to-last character: Hill function coefficient i.e. regulation "switchy-ness", choose from {2,3} or OMIT if =1
%          * last character: choose from {p,n} to specify positive 'p' regulation (promote parameter) or negative 'n' regulation (inhibit parameter)
%          Can specify multiple types of regulation by separating each type with underscore '_'
%          e.g. 'pG2p', 'pGp', 'vB1n_pG2p', 'vB22n_pG2p', 'none'
% OUTPUT:
%      funs: struct of function handles for each beta-cell parameter that
%         may be regulated. Each function handle is either a constant value or
%         a Hill function depending on whether a parameter is regulated.
%         This struct is input for functions such as ODE_RegulateBetaCells.m



% initialize function handles assuming no regulation of these parameters
nfun = @(G,I,B1,B2,par,KG,KI,KB)  par;
vfun = @(G,I,B1,B2,par,KG,KI,KB)  par;
pfun = @(G,I,B1,B2,par,KG,KI,KB)  par;
qfun = @(G,I,B1,B2,par,KG,KI,KB)  par;
mfun = @(G,I,B1,B2,par,KG,KI,KB)  par;

eachRegType = split(model_structure,'_'); % underscore separates each type of regulation
df = dict_fxn; dfp = dict_fxn_p; % get dictionaries of function handles

% input 'none' returns struct of function handles assuming no regulation
if strcmp(model_structure,'none')
    funs = struct('neo',nfun,'v',vfun,'p',pfun,'m',mfun,'q',qfun);
    return
end

for ms = 1:length(eachRegType) % for each regulation type, determine Hill coeff, get function handle
    regType = eachRegType{ms}; 
    fun = df(regType(2:end)); % get function handle from dictionary
    if strcmp(regType(1),'p') % if regulating p (self-renewal probability)...
        fun = dfp(regType(2:end)); % ... then multiply Hill (in [0,1]) by 1 instead of p0
    end
    switch regType(1) % switch depending on parameter (first character)
        case 'n' % neogenesis
            nfun = fun;
        case 'v' % cell division 
            vfun = fun;
        case 'p' % self-renewal probability
            pfun = fun;
        case 'q' % dedifferentiation
            qfun = fun; 
        case 'm' % mitosis-independent maturation
            mfun = fun;
    end
end

% output
funs = struct('neo',nfun,'v',vfun,'p',pfun,'m',mfun,'q',qfun);


end

%% ----- dictionaries for Hill function handles -----

function df = dict_fxn
    
    df = containers.Map();
    
    % regulation mediated by mature beta-cells (B2)
    df('B23n') = @(G,I,B1,B2,par,KG,KI,KB) par* KB^3 ./ (KB^3 + B2.^3);
    df('B23p') = @(G,I,B1,B2,par,KG,KI,KB) par* B2.^3 ./ (KB^3 + B2.^3);
    df('B22n') = @(G,I,B1,B2,par,KG,KI,KB) par* KB^2 ./ (KB^2 + B2.^2);
    df('B22p') = @(G,I,B1,B2,par,KG,KI,KB) par* B2.^2 ./ (KB^2 + B2.^2);
    df('B2n') = @(G,I,B1,B2,par,KG,KI,KB) par* KB ./ (KB + B2);
    df('B2p') = @(G,I,B1,B2,par,KG,KI,KB) par* B2 ./ (KB + B2);
    
    % regulation mediated by immature beta-cells (B1)
    df('B13n') = @(G,I,B1,B2,par,KG,KI,KB) par* KB^3 ./ (KB^3 + B1.^3);
    df('B13p') = @(G,I,B1,B2,par,KG,KI,KB) par* B1.^3 ./ (KB^3 + B1.^3);
    df('B12n') = @(G,I,B1,B2,par,KG,KI,KB) par* KB^2 ./ (KB^2 + B1.^2);
    df('B12p') = @(G,I,B1,B2,par,KG,KI,KB) par* B1.^2 ./ (KB^2 + B1.^2);
    df('B1n') = @(G,I,B1,B2,par,KG,KI,KB) par* KB ./ (KB + B1);
    df('B1p') = @(G,I,B1,B2,par,KG,KI,KB) par* B1 ./ (KB + B1);
    
    % regulation mediated by glucose
    df('G3n') = @(G,I,B1,B2,par,KG,KI,KB) par* KG^3 ./ (KG^3 + G.^3);
    df('G3p') = @(G,I,B1,B2,par,KG,KI,KB) par* G.^3 ./ (KG^3 + G.^3);
    df('G2n') = @(G,I,B1,B2,par,KG,KI,KB) par* KG^2 ./ (KG^2 + G.^2);
    df('G2p') = @(G,I,B1,B2,par,KG,KI,KB) par* G.^2 ./ (KG^2 + G.^2);
    df('Gn') = @(G,I,B1,B2,par,KG,KI,KB) par* KG ./ (KG + G);
    df('Gp') = @(G,I,B1,B2,par,KG,KI,KB) par* G ./ (KG + G);

    % regulation mediated by insulin
    df('I3n') = @(G,I,B1,B2,par,KG,KI,KB) par* KI^3 ./ (KI^3 + I.^3);
    df('I3p') = @(G,I,B1,B2,par,KG,KI,KB) par* I.^3 ./ (KI^3 + I.^3);
    df('I2n') = @(G,I,B1,B2,par,KG,KI,KB) par* KI^2 ./ (KI^2 + I.^2);
    df('I2p') = @(G,I,B1,B2,par,KG,KI,KB) par* I.^2 ./ (KI^2 + I.^2);
    df('In') = @(G,I,B1,B2,par,KG,KI,KB) par* KI ./ (KI + I);
    df('Ip') = @(G,I,B1,B2,par,KG,KI,KB) par* I ./ (KI + I);
    
    
end





function df = dict_fxn_p
    % Function handle dictionary for regulation of p (self-renewal probability) 
    % Hill function alone takes about the same value as p ~ 0.5 (and has same range in [0,1])
    % so Hill alone is better than multiplying p * Hill
    
    df = containers.Map();
   
    
    df('B23n') = @(G,I,B1,B2,par,KG,KI,KB)  KB^3 ./ (KB^3 + B2.^3);
    df('B23p') = @(G,I,B1,B2,par,KG,KI,KB)  B2.^3 ./ (KB^3 + B2.^3);
    df('B22n') = @(G,I,B1,B2,par,KG,KI,KB)  KB^2 ./ (KB^2 + B2.^2);
    df('B22p') = @(G,I,B1,B2,par,KG,KI,KB)  B2.^2 ./ (KB^2 + B2.^2);
    df('B2n') = @(G,I,B1,B2,par,KG,KI,KB)  KB ./ (KB + B2);
    df('B2p') = @(G,I,B1,B2,par,KG,KI,KB)  B2 ./ (KB + B2);
    
    df('B13n') = @(G,I,B1,B2,par,KG,KI,KB)  KB^3 ./ (KB^3 + B1.^3);
    df('B13p') = @(G,I,B1,B2,par,KG,KI,KB)  B1.^3 ./ (KB^3 + B1.^3);
    df('B12n') = @(G,I,B1,B2,par,KG,KI,KB)  KB^2 ./ (KB^2 + B1.^2);
    df('B12p') = @(G,I,B1,B2,par,KG,KI,KB)  B1.^2 ./ (KB^2 + B1.^2);
    df('B1n') = @(G,I,B1,B2,par,KG,KI,KB)  KB ./ (KB + B1);
    df('B1p') = @(G,I,B1,B2,par,KG,KI,KB)  B1 ./ (KB + B1);
    
    
    df('G3n') = @(G,I,B1,B2,par,KG,KI,KB)  KG^3 ./ (KG^3 + G.^3);
    df('G3p') = @(G,I,B1,B2,par,KG,KI,KB)  G.^3 ./ (KG^3 + G.^3);
    df('G2n') = @(G,I,B1,B2,par,KG,KI,KB)  KG^2 ./ (KG^2 + G.^2);
    df('G2p') = @(G,I,B1,B2,par,KG,KI,KB)  G.^2 ./ (KG^2 + G.^2);
    df('Gn') = @(G,I,B1,B2,par,KG,KI,KB)  KG ./ (KG + G);
    df('Gp') = @(G,I,B1,B2,par,KG,KI,KB)  G ./ (KG + G);

    df('I3n') = @(G,I,B1,B2,par,KG,KI,KB)  KI^3 ./ (KI^3 + I.^3);
    df('I3p') = @(G,I,B1,B2,par,KG,KI,KB)  I.^3 ./ (KI^3 + I.^3);
    df('I2n') = @(G,I,B1,B2,par,KG,KI,KB)  KI^2 ./ (KI^2 + I.^2);
    df('I2p') = @(G,I,B1,B2,par,KG,KI,KB)  I.^2 ./ (KI^2 + I.^2);
    df('In') = @(G,I,B1,B2,par,KG,KI,KB)  KI ./ (KI + I);
    df('Ip') = @(G,I,B1,B2,par,KG,KI,KB)  I ./ (KI + I);
    
    
end



