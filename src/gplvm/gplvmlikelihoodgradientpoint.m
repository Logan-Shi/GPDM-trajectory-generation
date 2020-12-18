function [L,g] = gplvmlikelihoodgradientpoint(x_i, i, A, invK, activeX, Y, D, theta, activeSet)

% GPLVMLIKELIHOODPOINT Compute gradient of data-point likelihood wrt x.

kbold = kernel(x_i, activeX, theta)';
f = A*kbold;
sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 

if sigma2 < 1e-6;
  sigma2 = 1e-6;
end

yHat_i = Y(i, :)' - f;

L = -D/2*log(sigma2) -D/2*log(2*pi) - yHat_i'*yHat_i/(2*sigma2)-x_i*x_i'/2;
L = -L;



if 0
% GPLVMGRADIENTPOINT Compute gradient of data-point likelihood wrt x.

kbold = kernel(x_i, activeX, theta)';
f = A*kbold;
sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 
yHat_i = Y(i, :)' - f;


if sigma2 < 1e-6;
  sigma2 = 1e-6;
end
end

prePart = theta(1)/sigma2*(yHat_i'*Y(activeSet, :)'+(D-yHat_i'*yHat_i/sigma2)*kbold')*invK;

g = zeros(size(x_i));
for k = 1:length(x_i)
  g(k) = prePart*((x_i(k) - activeX(:, k)).*kbold);
end

g = g + x_i;