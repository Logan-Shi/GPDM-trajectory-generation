function [X, theta, activeSet] = gplvmfit(X, Y, theta, numActive, optionsKernel, ...
			       options, extIters)

% GPLVMFIT Fit a Gaussian process latent variable model.

numData = size(Y, 1);
dataDim = size(Y, 2);

for iters = 1:extIters

  fprintf(2,'Iteration %d\n',iters);

  % Fit a gaussian proces to the model using ivm algorithm to sparsify.
  activeSet = gplvmivm(X, theta, numActive);

  % Optimise kernel parameters in log space
  lntheta = log(theta);

%  lntheta = fminunc('gplvmlikelihoodgradient', lntheta, optionsKernel, ...
%		Y(activeSet, :), X(activeSet, :));
  lntheta = scg('gplvmlikelihood', lntheta, optionsKernel, ...
		'gplvmgradient', Y(activeSet, :), X(activeSet, :));
  theta = exp(lntheta);

  fprintf('Kernel width: %4.2f\n', 1/theta(1))
  fprintf('RBF Process variance: %4.2f\n', theta(2))
  fprintf('Noise variance: %4.2f\n', 1/theta(3))
  
  % Update active set
  activeSet = gplvmivm(X, theta, numActive);
  
  % Compute new kernel and matrix of `alpha values'
  [K, invK] = computeKernel(X(activeSet, :), theta);
  A = Y(activeSet, :)'*invK;

  % Iterate through the data updating positions
  for i = 1:numData
    % do not update active set positions
    if ~any(activeSet==i)
%      X(i, :) = fminunc('gplvmlikelihoodgradientpoint', X(i, :),  options, ...
%		    i, A, invK, X(activeSet, :), ...
%		    Y, dataDim, theta, activeSet);
      X(i, :) = scg('gplvmlikelihoodpoint', X(i, :),  options, ...
		    'gplvmgradientpoint', i, A, invK, X(activeSet, :), ...
		    Y, dataDim, theta, activeSet);

    end
    if ~rem(i, floor(numData/10))
      fprintf('Finished point %d\n', i)
    end      
  end
end
