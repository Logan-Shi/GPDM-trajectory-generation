function [X, theta, thetap, weights] = gpdmfitFullX(X, Y, weights,...
    segments, theta, thetap, options, extIters, modelType, missing, fixX)

% learn full GPDM but fix hyperparameters

if (~exist('fixX', 'var'))
    fixX = 0;
end
N = size(Y,1);
D = size(Y,2);
q = size(X,2);

nmissing = setdiff(1:N, missing);

ndp = 0;
if (modelType(3) == 0)
	ndp = 4;
elseif (modelType(3) == 1)
	ndp = 3;
elseif (modelType(3) == 2)
	ndp = 2;
elseif (modelType(3) == 3)
    ndp = 6;
elseif (modelType(3) == 4)
    ndp = 3;
elseif (modelType(3) == 5)
    ndp = 4;
end

for iters = 1:extIters
    fprintf(2,'Iteration %d\n',iters);
    tic
    
    if (fixX == 1)
        lntheta = log(theta);
        lnthetap = log(thetap);
        
        params = [lntheta lnthetap];
        params = scg('gpdmlikelihood', params, options, 'gpdmgradientX',...
            Y, weights, segments, modelType, missing, X);
        
        lntheta = params(1:3);
        lnthetap = params(4:end);
        theta = exp(lntheta);
        thetap = exp(lnthetap);
    else
        lntheta = log(theta);
        lnthetap = log(thetap);
        
        params = [X(:)' lntheta lnthetap];
        params = scg('gpdmlikelihood', params, options, 'gpdmgradientX',...
            Y, weights, segments, modelType, missing);
        X = reshape(params(1:N*q), N, q);
        
        lntheta = params(end-(ndp+2):end-ndp);
        lnthetap = params(end-(ndp-1):end);
        theta = exp(lntheta);
        thetap = exp(lnthetap);
    end
    if (modelType(3) == 0)
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
        fprintf('Prior Kernel width1: %4.2f\n', 1/thetap(1))
        fprintf('Prior Kernel width2: %4.2f\n', 1/thetap(2))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(3))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(4))
    elseif (modelType(3) == 1)
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
        fprintf('Prior Kernel width: %4.2f\n', 1/thetap(1))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(2))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(3))
    elseif (modelType(3) == 2)
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
        fprintf('Prior Process variance: %4.2f\n', thetap(1))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(2))
    elseif (modelType(3) == 3)
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
        fprintf('Prior Linear Process variance1: %4.2f\n', thetap(1))
        fprintf('Prior Linear Process variance2: %4.2f\n', thetap(2))
        fprintf('Prior Kernel width1: %4.2f\n', 1/thetap(3))
        fprintf('Prior Kernel width2: %4.2f\n', 1/thetap(4))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(5))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(6))
    elseif (modelType(3) == 4)
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
        fprintf('Prior Process variance1: %4.2f\n', thetap(1))
        fprintf('Prior Process variance2: %4.2f\n', thetap(2))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(3))
    elseif (modelType(3) == 5)
        fprintf('Kernel width: %4.2f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.2f\n', theta(2))
        fprintf('Noise variance: %4.2f\n', 1/theta(3))
        fprintf('Prior Linear Process variance: %4.2f\n', thetap(1))
        fprintf('Prior Kernel width: %4.2f\n', 1/thetap(2))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(3))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(4))
    end
    
    if (fixX == 1)
        fprintf('-log Likelihood: %4.2f\n', ...
            gpdmlikelihood(params, Y, weights, segments, modelType, ...
			missing, X));
    else
        fprintf('-log Likelihood: %4.2f\n', ...
            gpdmlikelihood(params, Y, weights, segments, modelType, missing));
    end    
    
    [K, invK] = computeKernel(X(nmissing,:), theta);
    for d=1:D
        denom = Y(nmissing,d)'*invK*Y(nmissing,d);
        if (denom ~= 0.0)
            weights(d) = sqrt(length(nmissing)/denom);
        end
    end
    toc
end
																																					   


