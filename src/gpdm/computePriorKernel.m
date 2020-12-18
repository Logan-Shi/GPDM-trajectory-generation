function [K, invK, rbfK] = computePriorKernel(X, theta, type, Kn, invKn);

if nargin < 4
if (type == -1)
    K = eye(size(X,1), size(X,1));
else
    theta = thetaConstrain(theta);
    if (type == 0)
        K = kernel2(X, X, theta);
        rbfK = K; 
    elseif (type == 1)
        K = kernel(X, X, theta);
        rbfK = K; 
    elseif (type == 2)
        K = lin_kernel(X, X, theta);
        rbfK = eye(size(K)); 
    elseif (type == 3)
        rbfK = kernel2(X, X, theta(3:5)); 
        K = lin_kernel2(X, X, theta(1:2)) + rbfK;
    elseif (type == 4)
        K = lin_kernel2(X, X, theta);
    elseif (type == 5)
        rbfK = kernel(X, X, theta(2:3)); 
        K = lin_kernel(X, X, theta(1)) + rbfK;
    end
    K = K + eye(size(X, 1))*1/theta(end);
end
if nargout > 1
    invK = pdinv(K);
end
else
    n = size(Kn,1);
    Ainv = invKn;

    theta = thetaConstrain(theta);
    if (type == 0)
        B = kernel2(X(1:n,:), X(n+1:end,:), theta);
        rbfB = B; 
        rbfC = kernel2(X(n+1:end,:), X(n+1:end,:), theta); 
        C = rbfC + eye(size(X(n+1:end,:), 1))*1/theta(end);
        rbfKn = Kn - eye(n)*1/theta(end); 
    elseif (type == 1)
        B = kernel(X(1:n,:), X(n+1:end,:), theta);
        rbfB = B; 
        rbfC = kernel(X(n+1:end,:), X(n+1:end,:), theta); 
        C = rbfC + eye(size(X(n+1:end,:), 1))*1/theta(end);
        rbfKn = Kn - eye(n)*1/theta(end); 
    elseif (type == 2)
        B = lin_kernel(X(1:n,:), X(n+1:end,:), theta);
        C = lin_kernel(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);
        rbfB = zeros(size(B)); 
        rbfC = eye(size(C)); 
        rbfKn = eye(size(Kn)); 
    elseif (type == 3)
        rbfB = kernel2(X(1:n,:), X(n+1:end,:), theta(3:5)); 
        B = lin_kernel2(X(1:n,:), X(n+1:end,:), theta(1:2)) + rbfB;
        rbfC = kernel2(X(n+1:end,:), X(n+1:end,:), theta(3:5)); 
        C = lin_kernel2(X(n+1:end,:), X(n+1:end,:), theta(1:2)) + rbfC + ...
            eye(size(X(n+1:end,:), 1))*1/theta(end);
        rbfKn = Kn - lin_kernel2(X(1:n,:), X(1:n,:), theta(1:2)) - eye(n)*1/theta(end);
    elseif (type == 4)
        B = lin_kernel2(X(1:n,:), X(n+1:end,:), theta);
        C = lin_kernel2(X(n+1:end,:), X(n+1:end,:), theta) + eye(size(X(n+1:end,:), 1))*1/theta(end);
        rbfB = zeros(size(B)); 
        rbfC = eye(size(C)); 
        rbfKn = eye(size(Kn)); 
    elseif (type == 5)
        rbfB =  kernel(X(1:n,:), X(n+1:end,:), theta(2:3)); 
        B = lin_kernel(X(1:n,:), X(n+1:end,:), theta(1)) + rbfB;
        rbfC = kernel(X(n+1:end,:), X(n+1:end,:), theta(2:3)); 
        C = lin_kernel(X(n+1:end,:), X(n+1:end,:), theta(1)) + rbfC + ...
            eye(size(X(n+1:end,:), 1))*1/theta(end);
        rbfKn = Kn - lin_kernel(X(1:n,:), X(1:n,:), theta(1)) - eye(n)*1/theta(end); 
    end

    rbfK = [rbfKn rbfB; rbfB' rbfC]; 
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