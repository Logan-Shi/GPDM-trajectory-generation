function [K1] = lin_kernelDiffParams(arg1, arg2, arg3)

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

theta = thetaConstrain(theta);

K1 = [X1]*[X2]';
