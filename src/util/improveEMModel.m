function improveEMModel(input, output, extItr)

load(input);

opt = foptions;
opt(9) = 0;
opt(14) = 100;

nsamples = 100;
tau = 10;
hmcopt = foptions;
hmcopt(1) = 1;
hmcopt(7) = tau;
hmcopt(9) = 0;
hmcopt(14) = nsamples;
hmcopt(18) = 1/1500;

[X, theta, thetad, w, Xsamples] = gpdmfitEM(X, Y, w, segments, theta, thetad, ... 
        opt, hmcopt, extItr, modelType, output, Xsamples);

save(output, 'X', 'Y', 'w', 'theta', 'thetad', 'modelType', ...
'N', 'D', 'q', 'meanData', 'segments', 'initY', 'varY', 'missing', 'refY', 'Xinit', 'Xsamples');
