function [K1, K2, K3] = kernel2DiffParams(arg1, arg2, arg3);

% KERNELDIFFPARAMS Get gradients of kernel wrt its parameters.

if nargin < 3
  X1 = arg1;
  X2 = arg1;
  theta = arg2;
else
  X1 = arg1;
  X2 = arg2;
  theta = arg3;
end

q = size(X1,2)/2;

theta = thetaConstrain(theta);

K = kernel2(X1, X2, theta);
K3 = K/theta(3);
K2 = -0.5*dist2(X1(:,q+1:end), X2(:,q+1:end)).*K;
K1 = -0.5*dist2(X1(:,1:q), X2(:,1:q)).*K;
