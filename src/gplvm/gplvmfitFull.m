function [X, theta, activeSet, scale] = gplvmfitFull(X, Y, theta, numActive, optionsKernel, ...
			       options, extIters,scale)

% GPLVMFIT Fit a Gaussian process latent variable model.

numData = size(Y, 1);
dataDim = size(Y, 2);
latentDim = size(X, 2);

if ~exist('scale','var')
   scale = ones(dataDim,1);
   Yscaled = Y;
else
   Yscaled = zeros(size(Y));

  [K, invK] = computeKernel(X, theta);

    for d=1:dataDim
     scale(d) = 1/sqrt(Y(:,d)'*(invK/numData)*Y(:,d));
     Yscaled(:,d) = scale(d) * Y(:,d);
    end
end

for iters = 1:extIters

  fprintf(2,'Iteration %d\n',iters);

  % Optimise kernel parameters in log space
  lntheta = log(theta);

  params = [X(:)' lntheta];

  params = scg('gplvmlikelihood', params, optionsKernel, 'gplvmgradient', Yscaled);

  X = reshape(params(1:numData*latentDim), numData, latentDim);
  lntheta = params(numData*latentDim+1:end);
  theta = exp(lntheta);

  fprintf('Kernel width: %4.2f\n', 1/theta(1))
  fprintf('RBF Process variance: %4.2f\n', theta(2))
  fprintf('Noise variance: %4.2f\n', 1/theta(3))

  % Compute new kernel and matrix of `alpha values'
%  [K, invK] = computeKernel(X, theta);

  % update scale
  
%  for d=1:dataDim
%   scale(d) = 1/sqrt(Y(:,d)'*(invK/numData)*Y(:,d));
%   Yscaled(:,d) = scale(d) * Y(:,d);
%  end
end




activeSet = gplvmivm(X, theta, numActive);
