function kx = lin_kernel(x, xKern, theta)

% KERNEL Compute the linear kernel

kx = theta(1)*[x]*[xKern]';
