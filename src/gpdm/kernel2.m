function kx = kernel2(x, xKern, theta)

% KERNEL Compute the rbf kernel

q = size(x,2)/2;
n2 = dist2(x(:,1:q), xKern(:,1:q));
m2 = dist2(x(:,q+1:end), xKern(:,q+1:end));
rn = theta(1)/2;
rm = theta(2)/2;

kx = theta(3)*exp(-(n2*rn + m2*rm));
   