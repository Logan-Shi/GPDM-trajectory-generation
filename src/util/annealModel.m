function [theta, thetad] = annealModel(input, output, options, extIters, recVar, dynVar)

load(input);
    
Y = Y + sqrt(recVar)*randn(size(Y));

theta = [1 1 1];
thetad = ones(size(thetad));

if (recVar ~= 0)
    for iters = 1:extIters
        fprintf(2,'Iteration %d\n',iters);
        
        lntheta = log(theta);
        
        params = [lntheta];
        params = scg('gplvmlikelihood', params, options, 'gplvmgradient',...
            Y, X);
        
        lntheta = params;
        theta = exp(lntheta);
        
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
    end
end
if (dynVar ~= 0)
    [thetad] = gpdmAnnealDyn(X, segments, thetad, options, extIters, modelType, ...
        dynVar);
end

save(output, 'X', 'Y', 'w', 'theta', 'thetad', ...
'modelType', 'N', 'D', 'q','meanData', 'segments', 'initY', 'varY');

