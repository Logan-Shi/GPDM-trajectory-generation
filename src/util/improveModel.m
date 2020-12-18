function improveModel(input, output, extItr)

addpath netlab gplvm;

load(input);

opt = foptions;
opt(9) = 0;
opt(14) = 10;
if (~exist('missing', 'var'))
	missing = [];
end
[X theta thetad w] = ...
gpdmfitFull(X, Y, w, segments, theta, thetad, opt, extItr, modelType, ...
missing);
save(output, 'X', 'Y', 'w', 'theta', 'thetad', 'modelType', ...
'N', 'D', 'q', 'meanData', 'segments', 'initY', 'varY', 'missing', 'refY');
