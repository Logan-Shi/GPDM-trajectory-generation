function [A, B] = computeCondPriorKernel(X, theta, type, Kn);

n = size(Kn,1);

if (type == 0)
    A = kernel2(X(1:n,:), X(n+1:end,:), theta);
    B = kernel2(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);
elseif (type == 1)
    A = kernel(X(1:n,:), X(n+1:end,:), theta);
    B = kernel(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);
elseif (type == 2)
    A = lin_kernel(X(1:n,:), X(n+1:end,:), theta);
    B = lin_kernel(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);
elseif (type == 3)
    A = lin_kernel2(X(1:n,:), X(n+1:end,:), theta(1:2)) + kernel2(X(1:n,:), X(n+1:end,:), theta(3:5));
    B = lin_kernel2(X(n+1:end,:), X(n+1:end,:), theta(1:2)) + kernel2(X(n+1:end,:), X(n+1:end,:), theta(3:5)) + ...
        eye(size(X(n+1:end,:), 1))*1/theta(end);
elseif (type == 4)
    A = lin_kernel2(X(1:n,:), X(n+1:end,:), theta);
    B = lin_kernel2(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);
elseif (type == 5)
    A = lin_kernel(X(1:n,:), X(n+1:end,:), theta(1)) + kernel(X(1:n,:), X(n+1:end,:), theta(2:3));
    B = lin_kernel(X(n+1:end,:), X(n+1:end,:), theta(1)) + kernel(X(n+1:end,:), X(n+1:end,:), theta(2:3)) + ...
        eye(size(X(n+1:end,:), 1))*1/theta(end);
end
