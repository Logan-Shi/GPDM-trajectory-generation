function kx = lin_kernel2(x, xKern, theta)

% KERNEL Compute the linear kernel
q = size(x,2)/2;
kx = theta(1)*x(:,1:q)*xKern(:,1:q)' + theta(2)*x(:,q+1:end)*xKern(:,q+1:end)';