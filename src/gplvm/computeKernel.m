function [K, invK] = computeKernel(X, theta, Kn, invKn);

% COMPUTEKERNEL Compute the kernel matrix for data X with parameters theta.

if (nargin < 4)
    theta = thetaConstrain(theta);
    K = kernel(X, X, theta);
    K = K + eye(size(X, 1))*1/theta(end);
    n = size(K,1);
    if nargout > 1
        invK = pdinv(K);
    end
else
    n = size(Kn,1);
    Ainv = invKn;

    theta = thetaConstrain(theta);
    B = kernel(X(1:n,:), X(n+1:end,:), theta);
    C = kernel(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);

    K = [Kn B; B' C];


    if nargout > 1
        D = C - B' * Ainv * B;
        Dinv = pdinv(D);

        AA = Ainv + Ainv * B * Dinv * B' * Ainv;
        BB = - Ainv * B * Dinv;
        CC = Dinv;

        invK = [AA, BB; BB' CC];
    end
end
