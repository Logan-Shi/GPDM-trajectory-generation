function [A, B] = computeCondKernel(X, theta, Kn);

n = size(Kn,1);

A = kernel(X(1:n,:), X(n+1:end,:), theta);
B = kernel(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);


