function mstring = format_model_string(ms)
% Use in plot titles or labels with options 'Interpreter,'latex' 
% INPUT: string identifying type(s) of beta-cell regulation
%    e.g. 'pGp' 'pG2p' 'pG2p_vB1n' 'vIp' 
% NOTE: get_regulatory_functions.m will help you understand how strings will map to regulatory functions


if strcmp(ms,'none')
    mstring = ms; 
    return 
end


splt = split(ms,'_'); % underscore separates each type of regulation
mstring = ''; % initialize output string

for i = 1:length(splt) % for each regulation type
    ms_i = splt(i); ms_i = ms_i{1}; % get string for ith regulation type
    posneg = ms_i(end); % 'p' positive regulation (promote) or 'n' negative regulation (inhibit)
    param = ms_i(1); % parameter being regulated (e.g. 'p','v','q')

    if strcmp(ms_i(2),'B') 
        % if beta-cells (B1 or B2) mediate regulation
        dynamicVar = sprintf('\\beta_%s',  ms_i(3));  %dynamicVar = sprintf('%s_%s', ms_i(2), ms_i(3));
        disp(dynamicVar)
    else
        % glucose (G) or insulin (I) mediate regulation
        dynamicVar = ms_i(2); 
    end

    if strcmp(dynamicVar(end),ms_i(end-1))
        % if no exponent specified
        expo = 0;
    else
        expo = ms_i(end-1);
    end
    mstring_i = do_LaTeX_format(dynamicVar, expo, posneg, param);
    mstring = sprintf('%s(%s)',mstring,mstring_i); % append to mstring
end


end




%% ----- function to format string --------

function mt = do_LaTeX_format(dynamicVar, expo, posneg, param)
    
    % --- transform p/n into promote/inhibit symbols
    switch posneg
        case 'p'
            pn = '\rightarrow';
        case 'n'
            pn = '\dashv';
        otherwise
            pn = '?';
    end

    % --- exponentiate (Hill coeff > 1)?
    if expo
        mt = sprintf('$%s^%s %s %s$',dynamicVar, expo, pn, param);
    else
        mt = sprintf('$%s %s %s$',dynamicVar, pn, param);
    end

end






