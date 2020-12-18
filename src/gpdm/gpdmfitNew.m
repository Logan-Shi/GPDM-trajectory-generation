function [newX] = gpdmfitNew(newX, newY, missing, X, Y, w, ...
    segments, theta, thetap, options, extIters, modelType, slow_way)

if (~exist('slow_way', 'var'))
    slow_way = 0; 
end

M = size(newX,1);
D = size(Y,2);
q = size(X,2);

% [K invK] = computeKernel(X, theta);
% [Xin Xout] = priorIO(X, segments, modelType);
% [Kd invKd] = computePriorKernel(Xin, thetap, modelType(3));
% % gpdmconditional(newX(:)', [Y; newY], w, segments, modelType, ...
% %     missing, X, [theta thetap], K, invK, Kd, invKd)
% % gpdmconditionalgrad(newX(:)', [Y; newY], w, segments, modelType, ...
% %     missing, X, [theta thetap], K, invK, Kd, invKd)
% 
% options(9) = 1; 
% params = [newX(:)'];
% params = scg('gpdmconditional', params, options, 'gpdmconditionalgrad',...
%          [Y; newY], w, segments, modelType, missing, X, [theta thetap], K, invK, Kd, invKd);
% 
% %allX = [X; newX]; 
% %gpdmlikelihood([allX(:)' log(theta) log(thetap)], [Y; newY], w, [segments size(Y,1)+1], modelType) - ...
% %gpdmlikelihood([X(:)' log(theta) log(thetap)], Y, w, segments, modelType)
% 
% 
% 
% return; 

[K invK] = computeKernel(X, theta);
[Xin Xout] = priorIO(X, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetap, modelType(3));

for iters = 1:extIters
    fprintf(2,'Iteration %d\n',iters);
    
    params = [newX(:)'];
    if (slow_way == 1) 
         params = scg('gpdmlikelihood', params, options, 'gpdmgradient',...
         [Y; newY], w, segments, modelType, missing, X, [theta thetap]);
    else
        params = scg('gpdmlikelihood', params, options, 'gpdmgradient',...
         [Y; newY], w, segments, modelType, missing, X, [theta thetap], K, invK, Kd, invKd);
    end
    %gpdmlikelihood(params, [Y; newY], w, segments, modelType, missing, ...
	%X, [theta thetap])
    newX = reshape(params(1:M*q), M, q);
end
																																					   


