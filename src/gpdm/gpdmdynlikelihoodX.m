function Lp = gpdmdynlikelihoodX(params, fixedX, thetap, Km, invKm, segments, priorModel, initX)


% for optimizing new data only based on dynamics 

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
    else 
        X = [fixedX; initX(1,:); reshape(params(1:M*q), M, q)];
    end
else
    X = [fixedX; reshape(params(1:M*q), M, q)];
end

[Xin Xout] = priorIO(X, seg, priorModel);

Np = size(Xin, 1);
[Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3), Km, invKm);
%[Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));

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


