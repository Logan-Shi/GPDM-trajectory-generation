function g = gpdmdyngradientX(params, fixedX, thetap, Km, invKm, segments, priorModel, initX)

if (~exist('initX', 'var'))
    noInitX = 1;
else
    noInitX = 0;
end

noModel = 0;
% if we want to sample from the prior before any data
if (size(fixedX,1) == 1)
    noModel = 1;
end

ndp = 0;
if (priorModel(3) == 0)
	ndp = 4;
elseif (priorModel(3) == 1)
	ndp = 3;
elseif (priorModel(3) == 2)
	ndp = 2;
elseif (priorModel(3) == 3)
    ndp = 6;
elseif (priorModel(3) == 4)
    ndp = 3;
elseif (priorModel(3) == 5)
    ndp = 4;
end
q = size(fixedX, 2);
seg = [segments size(fixedX,1)+1];
if (noModel == 1)
    fixedX = [];
    seg = [1];
end
M = length(params)/q;


if (noInitX == 0)
    if (priorModel(1) < 2)
        X = [fixedX; initX(1,:); initX(2,:); reshape(params(1:M*q), M, q)];
        N = M + size(fixedX,1) + 2;
    else 
        X = [fixedX; initX(1,:); reshape(params(1:M*q), M, q)];
        N = M + size(fixedX,1) + 1;
    end
else
    X = [fixedX; reshape(params(1:M*q), M, q)];
    N = M + size(fixedX,1);
end

[Xin Xout] = priorIO(X, seg, priorModel);
Np = size(Xin, 1);
if (priorModel(1) == 0 || priorModel(1) == 1)
    Mp = M;
else
    Mp = M + 1;
end

qp = size(Xin, 2);

[Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3), Km, invKm);
%[Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));
dLp_dK = -q/2*invKp + .5*invKp*Xout*Xout'*invKp;

dLp_dxin = zeros(Np, qp);

for d = 1:qp
    if (priorModel(3) == 0)
        Kpart = kernel2DiffX(Xin, thetap, d);
    elseif (priorModel(3) == 1)
        Kpart = kernelDiffX(Xin, thetap, d);
    elseif (priorModel(3) == 2)
        Kpart = lin_kernelDiffX(Xin, thetap, d);
    elseif (priorModel(3) == 3)
        Kpart = lin_kernel2DiffX(Xin, thetap(1:2), d) + ...
            kernel2DiffX(Xin, thetap(3:5), d);
    elseif (priorModel(3) == 4)
        Kpart = lin_kernel2DiffX(Xin, thetap, d);
    elseif (priorModel(3) == 5)
        Kpart = lin_kernelDiffX(Xin, thetap(1), d) + ...
            kernelDiffX(Xin, thetap(2:3), d);
    end
    dLp_dxin(:, d) = ...
        2*sum(dLp_dK.*Kpart, 2) - diag(dLp_dK).*diag(Kpart);
end

dLp_dx = priorDiffX(dLp_dxin, -invKp*Xout, N, q, seg, priorModel);

if (priorModel(1) == 0 || priorModel(1) == 1)
    dLp_dx(seg+1,:) = dLp_dx(seg+1,:) - (X(seg+1,:) - X(seg,:));
    dLp_dx(seg,:) = dLp_dx(seg,:) - X(seg,:) + (X(seg+1,:) - X(seg,:));
elseif (priorModel(1) == 2 || priorModel(1) == 3)
    dLp_dx(seg,:) = dLp_dx(seg,:) - X(seg,:);
end

gX = dLp_dx(N-M+1:end,:);
g = -[gX(:)'];
