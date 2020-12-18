function KPart = lin_kernelDiffX(X, theta, q);

% KERNELDIFFX Compute the gradient of the kernel with respect to X.

% Since the result is all zeros apart from one row/column combination,
% and the result is symetric, here we return a matrix whose rows are the
% vectors which give the row/column combinations  

theta = thetaConstrain(theta);

numData = size(X, 1);

K = repmat(X(:, q)', numData, 1);

KPart = theta(1)*(tril(K) + triu(K));

