function [X, theta, thetap, weights, thetaSamples] = gpdmfitEM2(X, Y, weights,...
    segments, theta, thetap, options, hmcoptions, extIters, modelType, output, thetaSamples)

global MARGINAL_W
global W_VARIANCE

N = size(Y,1);
D = size(Y,2);
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

opt14 = options(14);
energyLog = []; 
done = 0; 

for iters = 1:extIters
    fprintf(2,'Iteration %d\n',iters);
    tic

    % E-step

    % Find an optimal X wrt current theta and thetad.
    options(14) = 10000;
    if (~exist('thetaSamples','var'))
        lntheta = log(theta);
        lnthetap = log(thetap);
        params = [lntheta lnthetap];
        
        params = scg('gpdmlikelihood', params, options, 'gpdmgradient',...
            Y, weights, segments, modelType, [], X);
    else
        params = mean(thetaSamples);
        
        params = scg('gpdmlikelihood', params, options, 'gpdmgradient',...
            Y, weights, segments, modelType, [], X);
    end
    theta = exp(params(1:3)); 
    thetap = exp(params(4:3+ndp)); 

    options(14) = opt14;

    initparams = params; 
    while (1)
        hmcoptions(14) = 1;
        hmcoptions(15) = 0; 
        % Sample from the posterior.
        params = initparams; 
        [thetaSamples, energies, diagn] = hmc('gpdmlikelihood', params, hmcoptions, 'gpdmgradient',...
            Y, weights, segments, modelType, [], X);
        if diagn.acc(1) < 1e-4
            ratio = 0.1;
        else
            hmcoptions(14) = 25;
            hmcoptions(15) = 10;
            % Sample from the posterior.
            params = initparams;
            [thetaSamples, energies, diagn] = hmc('gpdmlikelihood', params, hmcoptions, 'gpdmgradient',...
                Y, weights, segments, modelType, [], X);

            diagn.acc(find(diagn.acc > 1)) = 1;
            ratio = sum(diagn.acc)/size(diagn.acc, 1)
            hmcoptions(18)
        end
        if (ratio < 0.6)
            hmcoptions(18) = 0.7*hmcoptions(18);
        end
        if (ratio > 0.8)
            hmcoptions(18) = 1.5*hmcoptions(18);
        end
        if (ratio > 0.6 && ratio < 0.9)
            hmcoptions(15) = 0; 
            hmcoptions(14) = 25;

            [thetaSamples2] = hmc('gpdmlikelihood', thetaSamples(end,:), hmcoptions, 'gpdmgradient',...
                Y, weights, segments, modelType, [], X);
            thetaSamples = [thetaSamples; thetaSamples2];
            break;
        end
    end
    
    save(strcat(output, num2str(iters)), 'X', 'Y', 'weights', 'segments', 'theta', 'thetap', ...
        'thetaSamples', 'options', 'hmcoptions', 'modelType', 'energyLog');
thetaSamples = [log(theta) log(thetap)]; 
thetaSamples = [thetaSamples; thetaSamples]; 
    if (done == 1)
        return;
    end


    % M-step optimize theta, thetap wrt the incomplete likelihood.
    first = 1; 
    lastWeights = weights; 
    lastX = X;
    initWeights = weights;
    initX = lastX;
    changeW = 0;
    changeX = 0;
    while (1)
      if (MARGINAL_W ~= 1)
        R = size(thetaSamples,1);
        invK = zeros(N,N);
        for r = 1:R
            thetaR = thetaSamples(r,:);
            [Kr, invKr] = computeKernel(X, exp(thetaR(1:3)));
            invK = invK + invKr;
        end
        invK = invK/R;
      end
        lastF = options(8);
        if (first == 1 && iters > 1) 
            lastF
            params = [X(:)']; 
            thisF = gpdmincompletelikelihood2(params, Y, weights, segments, modelType, thetaSamples)
%             if (lastF - thisF < options(3))
%                 return;
%             end
        end
        if (MARGINAL_W ~= 1)
        for d=1:D
            denom = Y(:,d)'*invK*Y(:,d) + 1/W_VARIANCE;
            %if (denom ~= 0.0)
            weights(d) = sqrt(N/denom);
            %end
        end
        end
        changeW = max(abs(weights - lastWeights));
        lastWeights = weights;
        params = [X(:)']; 
        [params, options, flog, pointlog, scalelog] = scg('gpdmincompletelikelihood2', params, options, 'gpdmincompletelikelihoodgrad2',...
            Y, weights, segments, modelType, thetaSamples);        
        X = reshape(params(1:N*q), N, q); 
        changeX = max(sqrt(sum((X - lastX)'.^2)));
        lastX = X;
        
        if(size(flog,2) == options(14) && flog(end) == flog(1))
            options(14) = options(14) + 10;
        end
     
        if (first == 0 && abs(lastF - options(8)) < options(3) && ...
                changeW < options(2) && changeX < options(2) && ...
                size(flog,2) < options(14))

            break;
        end
        first = 0; 
    end
    changeW = max(abs(weights - initWeights));
    changeX = max(abs(X - initX));
    
    if (changeW < options(2) && changeX < options(2))
        done = 1; 
    end
    

    if (modelType(3) == 0)
        fprintf('Kernel width: %4.8f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.8f\n', theta(2))
        fprintf('Noise variance: %4.8f\n', 1/theta(3))
        fprintf('Prior Kernel width1: %4.8f\n', 1/thetap(1))
        fprintf('Prior Kernel width2: %4.8f\n', 1/thetap(2))
        fprintf('Prior RBF Process variance: %4.8f\n', thetap(3))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(4))
    elseif (modelType(3) == 1)
        fprintf('Kernel width: %4.8f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.8f\n', theta(2))
        fprintf('Noise variance: %4.8f\n', 1/theta(3))
        fprintf('Prior Kernel width: %4.8f\n', 1/thetap(1))
        fprintf('Prior RBF Process variance: %4.8f\n', thetap(2))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(3))
        fprintf('Prior signal to noise ratio: %4.8f\n', sqrt(thetap(2)*thetap(3))); 
        fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/thetap(1))); 
    elseif (modelType(3) == 2)
        fprintf('Kernel width: %4.8f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.8f\n', theta(2))
        fprintf('Noise variance: %4.8f\n', 1/theta(3))
        fprintf('Prior Process variance: %4.8f\n', thetap(1))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(2))
    elseif (modelType(3) == 3)
        fprintf('Kernel width: %4.8f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.8f\n', theta(2))
        fprintf('Noise variance: %4.8f\n', 1/theta(3))
        fprintf('Prior Linear Process variance1: %4.8f\n', thetap(1))
        fprintf('Prior Linear Process variance2: %4.8f\n', thetap(2))
        fprintf('Prior Kernel width1: %4.8f\n', 1/thetap(3))
        fprintf('Prior Kernel width2: %4.8f\n', 1/thetap(4))
        fprintf('Prior RBF Process variance: %4.8f\n', thetap(5))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(6))
    elseif (modelType(3) == 4)
        fprintf('Kernel width: %4.8f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.8f\n', theta(2))
        fprintf('Noise variance: %4.8f\n', 1/theta(3))
        fprintf('Prior Process variance1: %4.8f\n', thetap(1))
        fprintf('Prior Process variance2: %4.8f\n', thetap(2))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(3))
    elseif (modelType(3) == 5)
        fprintf('Kernel width: %4.8f\n', 1/theta(1))
        fprintf('RBF Process variance: %4.8f\n', theta(2))
        fprintf('Noise variance: %4.8f\n', 1/theta(3))
        fprintf('Prior Linear Process variance: %4.8f\n', thetap(1))
        fprintf('Prior Kernel width: %4.8f\n', 1/thetap(2))
        fprintf('Prior RBF Process variance: %4.8f\n', thetap(3))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(4))
    end

%     fprintf('-log incomplete Likelihood: %4.8f\n', ...
%         gpdmincompletelikelihood(params, Y, weights, segments, modelType, Xsamples, q));
    fprintf('-log incomplete Likelihood: %4.8f\n', options(8));

    energyLog = [energyLog options(8)]; 
    toc

%     if (size(energyLog, 2) > 1)
%         if ((energyLog(end-1) - energyLog(end)) < 1e-4)
%             return;
%         end
%     end
end


