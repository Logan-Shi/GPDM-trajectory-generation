function [X, theta, thetap, weights, Xsamples, I, Isegments] = gpdmfitEM(X, Y, weights,...
    segments, theta, thetap, options, hmcoptions, extIters, modelType, output, Xsamples)

global MARGINAL_W
global MARGINAL_DW
global W_VARIANCE
global SUBSET_SIZE
global FIX_HP

if MARGINAL_W == 1
    numInnerItr = 1;
else
    numInnerItr = 10;
end
I = 1:size(Y,1); 
Isegments = segments; 
fY = Y;
fX = X;
fsegments = segments;

if SUBSET_SIZE > 0
    numInnerItr = 10;
    [I Isegments] = selectSubset(Y, segments, SUBSET_SIZE); 

    Y = Y(I,:); 
    X = X(I,:); 
    segments = Isegments; 
    plotseries(X, segments, 'b'); 
    pause; 
end

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
num_samples = hmcoptions(14); 
energyLog = []; 
changeThetaLog = []; 
changeWLog = []; 
done = 0; 



for iters = 1:extIters
%     if (iters == extIters)
%         numInnerItr = 100; 
%     end
    fprintf(2,'Iteration %d\n',iters);
    tic
%     if SUBSET_SIZE > 0
%         [I Isegments] = selectSubset(fY, fsegments, SUBSET_SIZE);
% 
%         Y = fY(I,:);
%         X = fX(I,:);
%         segments = Isegments;
%         N = size(Y,1);
%     end

    % E-step

    % Find an optimal X wrt current theta and thetad.
    options(14) = 10000;
    if (hmcoptions(14) < 0)
        options(14) = 100; % We just do coordinate descent here.
    end
%     if SUBSET_SIZE > 0
%         params = [X(:)'];
%         params = scg('gpdmposterior', params, options, 'gpdmposteriorgrad',...
%             Y, weights, segments, modelType, [theta thetap]);
%     else
        if (~exist('Xsamples','var'))
            params = [X(:)'];
            params = scg('gpdmposterior', params, options, 'gpdmposteriorgrad',...
                Y, weights, segments, modelType, [theta thetap]);
        else
            if (size(Xsamples,1) == 1)
                params = Xsamples; 
            else
                params = mean(Xsamples);
            end
            params = scg('gpdmposterior', params, options, 'gpdmposteriorgrad',...
                Y, weights, segments, modelType, [theta thetap]);
        end
    %end
    Xsamples = params; 
    X = reshape(params(1:N*q), N, q);
    %     plotseries(X, segments, 'b');
    %     pause;
    options(14) = opt14;

    if (hmcoptions(14) > 0) 
    initparams = params; 
    while (1)
        hmcoptions(14) = 1;
        hmcoptions(15) = 0; 
        % Sample from the posterior.
        params = initparams; 
        [Xsamples, energies, diagn] = hmc('gpdmposterior', params, hmcoptions, 'gpdmposteriorgrad',...
            Y, weights, segments, modelType, [theta thetap]);
        if diagn.acc(1) < 1e-4
            ratio = 0.1;
        else
            hmcoptions(14) = 25;
            hmcoptions(15) = 10;
            % Sample from the posterior.
            params = initparams;
            [Xsamples, energies, diagn] = hmc('gpdmposterior', params, hmcoptions, 'gpdmposteriorgrad',...
                Y, weights, segments, modelType, [theta thetap]);

            diagn.acc(find(diagn.acc > 1)) = 1;
            ratio = sum(diagn.acc)/size(diagn.acc, 1)
            hmcoptions(18)
        end
        if (ratio < 0.6)
            hmcoptions(18) = 0.7*hmcoptions(18);
        end
        if (ratio > 0.9)
            hmcoptions(18) = 1.5*hmcoptions(18);
        end
        if (ratio > 0.6 && ratio < 0.95)
            hmcoptions(15) = 0; 
            hmcoptions(14) = num_samples - 25;
            if (hmcoptions(14) > 0) 
                [Xsamples2] = hmc('gpdmposterior', Xsamples(end,:), hmcoptions, 'gpdmposteriorgrad',...
                    Y, weights, segments, modelType, [theta thetap]);
                Xsamples = [Xsamples; Xsamples2];
            end
            break;
        end
    end
    hmcoptions(14) = num_samples; 
    
    save(strcat(output, num2str(iters)), 'X', 'Y', 'weights', 'segments', 'theta', 'thetap', ...
        'Xsamples', 'options', 'hmcoptions', 'modelType', 'energyLog', 'changeWLog', 'changeThetaLog');
    end
    if (done == 1)
        break;
    end


    % M-step optimize theta, thetap wrt the incomplete likelihood.
    first = 1; 
    lastWeights = weights; 
    lastTheta = [log(theta) log(thetap)]; 
    initWeights = weights;
    initTheta = lastTheta; 
    changeW = 0;
    changeTheta = 0;
    innerItr = 0; 
    while (1)
        innerItr = innerItr + 1;
        if (MARGINAL_W ~= 1)
            R = size(Xsamples,1);
            invK = zeros(N,N);
            for r = 1:R
                Xr = reshape(Xsamples(r,1:size(Xsamples,2)), size(Xsamples,2)/q, q);
                [Kr, invKr] = computeKernel(Xr, theta);
                invK = invK + invKr;
            end
            invK = invK/R;
        end
        lntheta = log(theta);
        lnthetap = log(thetap);
        lastF = options(8);
        params = [lntheta lnthetap];

        if (first == 1 && iters > 1)
            lastF
            thisF = gpdmincompletelikelihood(params, Y, weights, segments, modelType, Xsamples, q)
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
%         theta = thetaConstrain(theta);
%         thetap = thetaConstrain(thetap);
        [params, options, flog, pointlog, scalelog] = scg('gpdmincompletelikelihood', params, options, 'gpdmincompletelikelihoodgrad',...
            Y, weights, segments, modelType, Xsamples, q);
        lntheta = params(1:3);
        lnthetap = params(4:end);
        theta = exp(lntheta);
        thetap = exp(lnthetap);
%         theta = thetaConstrain(theta);
%         thetap = thetaConstrain(thetap);
        changeTheta = max(abs([log(theta) log(thetap)] - lastTheta));
        lastTheta = [log(theta) log(thetap)];

        if(size(flog,2) == options(14) && flog(end) == flog(1))
            options(14) = options(14) + 10;
        end

        if (first == 0 && abs(lastF - options(8)) < options(3) && ...
                changeW < options(2) && changeTheta < options(2) && ...
                size(flog,2) < options(14) || innerItr == numInnerItr)

            break;
        end
        first = 0;
    end
    changeW = max(abs(weights - initWeights));
    changeTheta = max(abs([log(theta) log(thetap)] - initTheta));

    if (changeW < options(2) && changeTheta < options(2))
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
    changeWLog = [changeWLog changeW]; 
    changeThetaLog = [changeThetaLog changeTheta]; 
    toc

    %     if (size(energyLog, 2) > 1)
    %         if ((energyLog(end-1) - energyLog(end)) < 1e-4)
    %             return;
    %         end
    %     end
end
N = size(fY,1);
if MARGINAL_W == 1
    R = size(Xsamples,1);
    invK = zeros(N,N);
    for r = 1:R
        Xr = reshape(Xsamples(r,1:size(Xsamples,2)), size(Xsamples,2)/q, q);
        [Xin Xout] = priorIO(Xr, segments, modelType);
        [Krp, invKrp] = computePriorKernel(Xin, thetap, modelType(3));
        [Kr, invKr] = computeKernel(Xr, theta);
        invK = invK + invKr;
        if r == 1
            invKp = invKrp;
        else
            invKp = invKp + invKrp;
        end
    end
    invK = invK/R;
    invKp = invKp/R;

    for d=1:D
        denom = Y(:,d)'*invK*Y(:,d) + 1/W_VARIANCE;
        weights(d) = sqrt(N/denom);
    end
end

if MARGINAL_DW == 1
    denom = 0;
    for d = 1:q
        denom = denom + Xout(:,d)'*invKp*Xout(:,d);
    end
    denom = denom + 1/W_VARIANCE;
    priorw = size(Xin,1)*q/denom;
    thetap(2) = thetap(2)/priorw;
    thetap(3) = thetap(3)*priorw;
end

if (modelType(3) == 1)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Kernel width: %4.8f\n', 1/thetap(1))
    fprintf('Prior RBF Process variance: %4.8f\n', thetap(2))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetap(3))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt(thetap(2)*thetap(3)));
    fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/thetap(1)));
end

if SUBSET_SIZE > 0
    FIX_HP = 1; 
    options(14) = 1000; 
    [X theta thetap weights] = gpdmfitFull(fX, fY, weights, fsegments, theta, thetap, options, ... 
     1, modelType, [], 0);
end

