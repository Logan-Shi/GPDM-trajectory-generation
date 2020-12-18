function g = gpdmdyngradientTheta(params, X, seg, priorModel, Xin, Xout)

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

q = size(X, 2);
Np = size(Xin, 1);
qp = size(Xin, 2);
lnthetap = params(end-(ndp-1):end);
thetap = exp(lnthetap);
thetap = thetaConstrain(thetap);


gParamp = zeros(1,ndp);


[Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));
dLp_dK = -q/2*invKp + .5*invKp*Xout*Xout'*invKp;

dk = zeros(1, ndp);

if (priorModel(3) == 0)
	[dK{1}, dK{2}, dK{3}] = kernel2DiffParams(Xin, thetap);
elseif (priorModel(3) == 1)
	[dK{1}, dK{2}] = kernelDiffParams(Xin, thetap);
elseif (priorModel(3) == 2)
	[dK{1}] = lin_kernelDiffParams(Xin, thetap);
elseif (priorModel(3) == 3)
	[dK{1} dK{2}] = lin_kernel2DiffParams(Xin, thetap(1:2));
	[dK{3}, dK{4}, dK{5}] = kernel2DiffParams(Xin, thetap(3:5));
elseif (priorModel(3) == 4)
	[dK{1} dK{2}] = lin_kernel2DiffParams(Xin, thetap);
elseif (priorModel(3) == 5)
	[dK{1}] = lin_kernelDiffParams(Xin, thetap(1));
	[dK{2}, dK{3}] = kernelDiffParams(Xin, thetap(2:3));
	end

for i = 1:(ndp-1)
	dk(i) = sum(sum(dLp_dK.*dK{i}));
end
dk(ndp) = -sum(diag(dLp_dK)/thetap(end).^2);

gParamp = dk.*thetap-1;

g = -gParamp;

