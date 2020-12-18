function L = gpdmconditional(params, Y, w, segments, priorModel, ...
    missing, Xn, fixedTheta, Kn, invKn, Kdn, invKdn)

global W_VARIANCE

N = size(Y, 1);
D = size(Y, 2);

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

q = size(Xn, 2);
M = N - size(Xn, 1);

Ym = Y(N-M+1:end, :);
Yn = Y(1:N-M, :); 
Xm = reshape(params(1:M*q), M, q); 
X = [Xn; Xm];
seg = [segments N-M+1];

[Xnin Xnout] = priorIO(Xn, segments, priorModel); 
[Xin Xout] = priorIO(X, seg, priorModel);
[Xmin Xmout] = priorIO(Xm, [1], priorModel);
nmissing = setdiff(1:M, missing);

Np = size(Xin, 1);
Mp = Np - size(Xnin, 1); 

theta = fixedTheta(end-(ndp+2):end-ndp);
thetap = fixedTheta(end-(ndp-1):end);

[A B] = computeCondKernel([Xn; Xm(nmissing,:)], theta, Kn);
Km = B - A'*invKn*A; 
invKm = pdinv(Km); 

numData = length(nmissing); 
L = -D*numData/2*log(2*pi)-D/2*logdet(Km);
for d= 1:D
    diff = Ym(nmissing, d) - A'*invKn*Yn(:,d); 
    L = L - 0.5*w(d)*w(d)*diff'*invKm*diff; 
end

L = L + numData*sum(log(w)); 

if (priorModel(3) == -1)
    Lp = -M*q/2*log(2*pi) -0.5*sum(sum(Xm.*Xm));
else
    [Ad Bd] = computeCondPriorKernel(Xin, thetap, priorModel(3), Kdn);
    Kdm = Bd - Ad'*invKdn*Ad; 
    invKdm = pdinv(Kdm); 
        Lp = -q*Mp/2*log(2*pi)-q/2*logdet(Kdm);

        for d= 1:q
            diff = Xmout(:, d) - Ad'*invKdn*Xnout(:,d);  
            Lp = Lp - 0.5*diff'*invKdm*diff;
        end

        if (priorModel(1) == 0 || priorModel(1) == 1)
            XDiff = Xm(2,:) - Xm(1,:);
            Lp = Lp - q/2*log(2*pi) - 0.5*sum(sum(XDiff.*XDiff));
        end

        Lp = Lp - q/2*log(2*pi) - 0.5*sum(sum(Xm(1,:).*Xm(1,:))); 
end

%L = L + Lp; 

L = -L;


