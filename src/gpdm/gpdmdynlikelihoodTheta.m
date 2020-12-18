function Lp = gpdmdynlikelihoodTheta(params, X, seg, priorModel, Xin, Xout)

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
[Xin Xout] = priorIO(X, seg, priorModel);
Np = size(Xin, 1);
lnthetap = params(end-(ndp-1):end);
thetap = exp(lnthetap);
thetap = thetaConstrain(thetap);

[Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));

Lp = -q*Np/2*log(2*pi)-q/2*logdet(Kp);

for d= 1:q
Lp = Lp -0.5*Xout(:, d)'*invKp*Xout(:, d);
end

if (priorModel(1) == 0 || priorModel(1) == 1)
	XDiff = X(seg+1,:) - X(seg,:);
	Lp = Lp - 0.5*sum(sum(XDiff.*XDiff));
end

Lp = Lp - 0.5*sum(sum(X(seg,:).*X(seg,:)));
Lp = Lp - sum(log(thetap));

Lp = -Lp;


