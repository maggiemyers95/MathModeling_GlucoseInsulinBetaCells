function rgb = v2rgb( v, vs )
% VECTOR TO RGB 
% Taken from http://stackoverflow.com/questions/11642826/use-matlab-colour-scheme-to-convert-float-to-rgb
% Inputs: v is an element in the array vs 

% normalize
minV = min(vs);
maxV = max(vs);
v = (v - minV)/(maxV - minV); 

% choose colormap e.g. 'jet' 'parula' 'jet' 'cool' 'spring' 'winter' 'copper' 'hsv'
cm = colormap('jet');

colorID = max(1, sum(v > [0:1/length(cm(:,1)):1])); % color ID has max value of 64

rgb = cm(colorID,:);

end
