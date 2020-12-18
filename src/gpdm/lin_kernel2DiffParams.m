function [K1, K2] = lin_kernel2DiffParams(arg1, arg2, arg3)

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

K1 = X1(:,1:q)*X2(:,1:q)';
K2 = X1(:,q+1:end)*X2(:,q+1:end)';
