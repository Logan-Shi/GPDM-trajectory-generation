function [X, theta, thetap, weights] = gpdmfitFull(X, Y, weights,...
    segments, theta, thetap, options, extIters, modelType, missing, fixX)

%global LOGTHETA
global USE_GAMMA_PRIOR
global GAMMA_ALPHA
global FIX_HP
global MARGINAL_W
global LEARN_SCALE
global MARGINAL_DW
global W_VARIANCE
global USE_OLD_MISSING_DATA

if MARGINAL_W == 1
    weights = zeros(size(weights)); 
end


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

opt14 = options(14); 
lastWeights = weights; 
lastX = X; 
%if ~LOGTHETA 
    lastTheta = [log(theta) log(thetap)]; 
%else
%    lastTheta = [theta(1:end-1) log(theta(end)) thetap(1:end-1) log(thetap(end))]; 
%end
energyLog = []; 
changeW = 0;
changeX = 0; 
changeTheta = 0;
for iters = 1:extIters
    fprintf(2,'Iteration %d\n',iters);
  
    if (FIX_HP ~= 1) && (MARGINAL_W ~= 1) && (LEARN_SCALE == 1)
        [K, invK] = computeKernel(X(nmissing,:), theta);
        if (USE_OLD_MISSING_DATA == 1) 
        else
        if (~isempty(missing))
            kbold = kernel(X(missing,:), X(nmissing,:), theta)';
            A = Y(nmissing,:)'*invK*kbold;
            Y(missing,:) = A';
            nmissing = 1:N;
            [K, invK] = computeKernel(X(nmissing,:), theta);
        end
        end

        for d=1:D
            if (W_VARIANCE == 0)
                denom = Y(nmissing,d)'*invK*Y(nmissing,d); 
                if (denom == 0) 
                    weights(d) = 1; 
                else 
                    weights(d) = sqrt(length(nmissing)/denom); 
                end
            else

                denom = Y(nmissing,d)'*invK*Y(nmissing,d) + 1/W_VARIANCE;

                weights(d) = sqrt(length(nmissing)/denom);
                %             weights(d) = 1;
            end

        end
    end
        changeW = max(abs(weights - lastWeights)); 
    lastWeights = weights; 
%     theta = thetaConstrain(theta);
%     thetap = thetaConstrain(thetap);
    if (fixX == 1)
%         if ~LOGTHETA
            lntheta = log(theta);
            lnthetap = log(thetap);
%         else
%             lntheta = [theta(1:end-1) log(theta(end))];
%             lnthetap = [thetap(1:end-1) log(thetap(end))];
%         end
        params = [lntheta lnthetap];
        [params options flog] = scg('gpdmlikelihood', params, options, 'gpdmgradient',...
            Y, weights, segments, modelType, missing, X);
        
        lntheta = params(1:3);
        lnthetap = params(4:end);
%         if ~LOGTHETA
            theta = exp(lntheta);
            thetap = exp(lnthetap);
%         else
%             theta = [lntheta(1:end-1) exp(lntheta(end))];
%             thetap = [lnthetap(1:end-1) exp(lnthetap(end))];
%         end
    else
%         if ~LOGTHETA
            lntheta = log(theta);
            lnthetap = log(thetap);
%         else
%             lntheta = [theta(1:end-1) log(theta(end))];
%             lnthetap = [thetap(1:end-1) log(thetap(end))];
%         end
%         
        params = [X(:)' lntheta lnthetap];
        [params options flog] = scg('gpdmlikelihood', params, options, 'gpdmgradient',...
            Y, weights, segments, modelType, missing);
        
        X = reshape(params(1:N*q), N, q); 
        changeX = max(max(abs(X - lastX))); 
        lastX = X; 
        
        lntheta = params(end-(ndp+2):end-ndp);
        lnthetap = params(end-(ndp-1):end);      
        
%         if ~LOGTHETA
            theta = exp(lntheta);
            thetap = exp(lnthetap);
%         else
%             theta = [lntheta(1:end-1) exp(lntheta(end))];
%             thetap = [lnthetap(1:end-1) exp(lnthetap(end))];
%         end
  
    end
%     theta = thetaConstrain(theta);
%     thetap = thetaConstrain(thetap);

% 
%     if ~LOGTHETA
        changeTheta = max(abs([log(theta) log(thetap)] - lastTheta));
        lastTheta = [log(theta) log(thetap)];
%     else
%         changeTheta = max(abs([theta(1:end-1) log(theta(end)) ...
%             thetap(1:end-1) log(thetap(end))] - lastTheta));
%         lastTheta = [theta(1:end-1) log(theta(end)) thetap(1:end-1) log(thetap(end))];
%     end

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
        fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3))); 
        fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1))); 
        fprintf('Prior Kernel width: %4.8f\n', 1/thetap(1))
        fprintf('Prior RBF Process variance: %4.8f\n', thetap(2))
        fprintf('Prior Noise variance: %4.8f\n', 1/thetap(3))
        fprintf('Prior signal to noise ratio: %4.8f\n', sqrt(thetap(2)*thetap(3))); 
        fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/thetap(1))); 
%         [Xin Xout] = priorIO(X, segments, modelType);
%         [Kp] = computePriorKernel(Xin, thetap, modelType(3));
%         q/2*logdet(Kp)
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
    fprintf('Delta X: %4.8f\n', changeX); 
    fprintf('Delta W: %4.8f\n', changeW); 
    fprintf('Delta Theta: %4.8f\n', changeTheta); 
%     if (fixX == 1)
%         fprintf('-log Likelihood: %4.8f\n', ...
%             gpdmlikelihood(params, Y, weights, segments, modelType, ...
% 			missing, X));
%     else
%         L = gpdmlikelihood(params, Y, weights, segments, modelType, missing); 
%         fprintf('-log Likelihood: %4.8f\n', L);
%         table = [table L]; 
% 
%         
%     end    
    
    fprintf('-log Posterior: %4.8f\n', options(8));
%     fprintf('-log Likelihood: %4.8f\n', options(8) + sum(log(theta)) + ... 
%         sum(log(thetap)) + 0.5*sum(weights.*weights));
%     fprintf('-log Likelihood: %4.8f\n', options(8) + sum(log(theta)) + ... 
%         0.5*sum(weights.*weights) + ... 
%         0.5*(thetap(1) - 0.2)*(thetap(1) - 0.2) + 0.5*log(2*pi) + ...
%         0.5*(thetap(2) - 10000/thetap(3))*(thetap(2) - 10000/thetap(3)) + 0.5*log(2*pi) + ...
%         log(thetap(3)));
% if gammaPrior == 1
%                 phi = 1./(thetap.*thetap); 
%         phi(end) = 1/phi(end); 
%          omega(2) = thetap(3)*thetap(3)/100^4; 
%         
%     fprintf('-log Likelihood: %4.8f\n', options(8) + sum(log(theta)) + ... 
%          0.5*sum(weights.*weights) + sum(0.5*alpha.*log(0.5*alpha./omega) + ...
%             (0.5*alpha - 1).*log(phi) - log(gamma(0.5*alpha)) - ...
%             0.5*(phi.*alpha)./omega));  
% else
%          fprintf('-log Likelihood: %4.8f\n', options(8) + sum(log(theta)) + ... 
%          sum(log(thetap)) + 0.5*sum(weights.*weights));
% end

    energyLog = [energyLog options(8)];
    
    w = weights;
    thetad = thetap; 
    if (mod(iters,100) == 0) 
    save(strcat('lastM', num2str(iters)), 'X', 'Y', 'w', 'theta', ...
        'thetad', 'modelType', 'N', 'D', 'q', 'segments', 'missing', ...
        'options', 'energyLog');
    end

    if (size(energyLog, 2) > 1)
        if ((energyLog(end-1) - energyLog(end)) < options(3) && ...
                changeX < options(2) && changeW < options(2) && ...
                changeTheta < options(2) && size(flog,2) < options(14))
            changeX
            changeW
            changeTheta
            size(flog,2)
            options(14)
            
            return;
        end
    end
  
    if(size(flog,2) == options(14) && flog(end) == flog(1))
        options(14) = opt14 + 10;
    else
        options(14) = opt14; 
    end
    
    

%     plot(weights);
%     pause;
end
nmissing = 1:N; 
if MARGINAL_W == 1
    [K, invK] = computeKernel(X(nmissing,:), theta);
    for d=1:D
        denom = Y(nmissing,d)'*invK*Y(nmissing,d) + 1/W_VARIANCE;
        weights(d) = sqrt(length(nmissing)/denom);
    end
end

if modelType(3) == 1
if MARGINAL_DW == 1
    [Xin Xout] = priorIO(X, segments, modelType);
    [Kp, invKp] = computePriorKernel(Xin, thetad, modelType(3));
    denom = 0; 
    for d = 1:q
        denom = denom + Xout(:,d)'*invKp*Xout(:,d); 
    end
    denom = denom + 1/W_VARIANCE;
    priorw = size(Xin,1)*q/denom; 
    thetap(2) = thetap(2)/priorw;
    thetap(3) = thetap(3)*priorw; 
end
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


