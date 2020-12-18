function [thetap] = gpdmAnnealDyn(X, segments, thetap, options, ...
extIters, modelType, noise_var)

% only fit dynamic hyperparameters

N = size(X,1);
q = size(X,2);

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

[Xin Xout] = priorIO(X, segments, modelType);
Xout = Xout + sqrt(noise_var)*randn(size(Xout));

for iters = 1:extIters
    fprintf(2,'Iteration %d\n',iters);
    lnthetap = log(thetap);
    
    params = [lnthetap];
    params = scg('gpdmdynlikelihoodTheta', params, options, ...
	'gpdmdyngradientTheta', X, segments, modelType, Xin, Xout);
    lnthetap = params(end-(ndp-1):end);
    thetap = exp(lnthetap);
    
	if (modelType(3) == 0)
        fprintf('Prior Kernel width1: %4.2f\n', 1/thetap(1))
        fprintf('Prior Kernel width2: %4.2f\n', 1/thetap(2))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(3))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(4))
    elseif (modelType(3) == 1)
        fprintf('Prior Kernel width: %4.2f\n', 1/thetap(1))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(2))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(3))
    elseif (modelType(3) == 2)
        fprintf('Prior Process variance: %4.2f\n', thetap(1))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(2))
    elseif (modelType(3) == 3)
        fprintf('Prior Linear Process variance1: %4.2f\n', thetap(1))
        fprintf('Prior Linear Process variance2: %4.2f\n', thetap(2))
        fprintf('Prior Kernel width1: %4.2f\n', 1/thetap(3))
        fprintf('Prior Kernel width2: %4.2f\n', 1/thetap(4))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(5))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(6))
    elseif (modelType(3) == 4)
        fprintf('Prior Process variance1: %4.2f\n', thetap(1))
        fprintf('Prior Process variance2: %4.2f\n', thetap(2))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(3))
    elseif (modelType(3) == 5)
        fprintf('Prior Linear Process variance: %4.2f\n', thetap(1))
        fprintf('Prior Kernel width: %4.2f\n', 1/thetap(2))
        fprintf('Prior RBF Process variance: %4.2f\n', thetap(3))
        fprintf('Prior Noise variance: %4.2f\n', 1/thetap(4))
    end
end
																																					   


