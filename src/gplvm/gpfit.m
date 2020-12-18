function [theta, activeSet] = gpfit(X, Y, theta, numActive, optionsKernel, ...
			       options)  
numData = size(Y, 1);
dataDim = size(Y, 2);
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