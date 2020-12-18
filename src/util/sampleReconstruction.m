function [Ysample Ymean] = sampleReconstruction(X_star, X, Y, theta, w, Kn, invKn) 

N = size(X,1); 
M = size(X_star,1); 
D = size(Y,2); 
[A B] = computeCondKernel([X; X_star], theta, Kn);
Km = B - A'*invKn*A; 
mu = []; 
C = []; 
Ymean = zeros(M,1); 
rootKm = Km^(1/2); 

z = randn(M,1);
for q = 1:D
    Ymean(:,q) = A'*invKn*Y(:,q);
    Ysample(:,q) = (1/w(q))*rootKm*z + Ymean(:,q);
end



