function L = gpdmposterior(params, Y, w, segments, priorModel, fixedTheta)

global INFO
global MARGINAL_W
global MARGINAL_DW
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

q = length(params)/N;
X = reshape(params(1:N*q), N, q);
seg = segments;
[Xin Xout] = priorIO(X, seg, priorModel);

Np = size(Xin, 1);

theta = fixedTheta(end-(ndp+2):end-ndp);
thetap = fixedTheta(end-(ndp-1):end);
theta = thetaConstrain(theta);
thetap = thetaConstrain(thetap);

[K, invK] = computeKernel(X, theta);

numData = N; 
if MARGINAL_W
    L = D*log(2)-D*log(sqrt(W_VARIANCE))-(0.5*numData+0.5)*D*log(2*pi)-D/2*logdet(K);
    for d = 1:D
        varTerm = 0.5*Y(:, d)'*invK*Y(:, d)+0.5/W_VARIANCE;
        if mod(numData,2) == 0
            L = L + (numData/2)*(log(numData-1) - log(2*varTerm)) + 0.5*(log(pi) - log(4*varTerm));
        else
            L = L + floor(numData/2)*(log(numData-1) - log(2*varTerm)) - log(2*varTerm);
        end
    end
    L = L - sum(log(theta));
else
    L = -D*N/2*log(2*pi)-D/2*logdet(K) + N*sum(log(w)); 
    for d= 1:D
        L = L -0.5*w(d)*w(d)*Y(:, d)'*invK*Y(:, d);
    end
    if (INFO == 1)
        fprintf('-L reconstruction: %4.2f\n', -L);
    end
    L = L - sum(log(theta))  + D*log(2) - D/2*log(2*pi*W_VARIANCE) - 0.5/W_VARIANCE*sum(w.*w);
end

if (priorModel(3) == -1)
    Lp = -N*q/2*log(2*pi)-0.5*sum(sum(X.*X));
else
    [Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));
    if MARGINAL_DW == 1
        Nq = Np*q;
        Lp = log(2/sqrt(W_VARIANCE))-(0.5*N*q+0.5)*log(2*pi)-q/2*logdet(Kp);
        varTerm = 0;
        for d = 1:q
            varTerm = varTerm + 0.5*Xout(:, d)'*invKp*Xout(:, d);
        end
        varTerm = varTerm + 0.5/W_VARIANCE;
        if mod(Nq,2) == 0
            Lp = Lp + (Nq/2)*(log(Nq-1) - log(2*varTerm)) + 0.5*(log(pi) - log(4*varTerm));
        else
            Lp = Lp + floor(Nq/2)*(log(Nq-1) - log(2*varTerm)) - log(2*varTerm);
        end
        Lp = Lp - sum(log(thetap));
        Lp = Lp - size(seg,2)*q/2*log(2*pi) - 0.5*sum(sum(X(seg,:).*X(seg,:)));
    else

        Lp = -q*Np/2*log(2*pi)-q/2*logdet(Kp);

        if (INFO == 1)
            fprintf('-L dynamics normalization: %4.2f\n', -Lp);
        end

        for d= 1:q
            Lp = Lp - 0.5*Xout(:, d)'*invKp*Xout(:, d);
        end

        if (priorModel(1) == 0 || priorModel(1) == 1)
            XDiff = X(seg+1,:) - X(seg,:);
            Lp = Lp - size(seg,2)*q/2*log(2*pi) - 0.5*sum(sum(XDiff.*XDiff));
        end

        Lp = Lp - size(seg,2)*q/2*log(2*pi) - 0.5*sum(sum(X(seg,:).*X(seg,:)));

        Lp = Lp - sum(log(thetap));
    end

end

L = L + Lp;


if (INFO == 1)
    fprintf('-L dynamics: %4.2f\n', -Lp);
end


L = -L;
if (INFO == 1)
	fprintf('-Log Likelihood: %4.2f\n', L);
end

