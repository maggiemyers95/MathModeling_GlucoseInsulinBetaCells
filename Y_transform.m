function Y_out = Y_transform(y)
% model-specific

if size(y,1)==4 % if model includes beta-cell populations... (ODE_RegulateBetaCells.m)
    y(5,:) = y(3,:)+y(4,:); % ...create new row for total beta-cells...
    y(6,:) = 100.*y(4,:)./y(5,:); % ... and new row for percent B2 cells
else % e.g. for model with only glucose and insulin (ODE_GlucoseInsulin.m)
    numVars = size(y,1); 
    y((numVars+1):6,:) = zeros((6-numVars),size(y,2)); % fill rows with zeros to make matrix size compatible with other model  
end
Y_out = y;
end