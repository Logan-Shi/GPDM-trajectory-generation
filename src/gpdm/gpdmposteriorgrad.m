function g = gpdmposteriorgrad(params, Y, w, segments, priorModel, fixedTheta)

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
M = N;
seg = segments;
[Xin Xout] = priorIO(X, seg, priorModel);

Np = size(Xin, 1);
qp = size(Xin, 2);


theta = fixedTheta(end-(ndp+2):end-ndp);
thetap = fixedTheta(end-(ndp-1):end);
theta = thetaConstrain(theta);
thetap = thetaConstrain(thetap);

[K, invK] = computeKernel(X, theta);
cacheK = K - eye(size(X, 1))*1/theta(end);

numData = N; 
if MARGINAL_W
    dL_dK = -D/2*invK;
    for d = 1:D
        varTerm = 0.5*Y(:, d)'*invK*Y(:, d)+0.5/W_VARIANCE;
        dvarTerm = -0.5*invK*Y(:, d)*Y(:, d)'*invK;
        if mod(numData,2) == 0
            dL_dK = dL_dK - ((numData/2)/varTerm + 0.5/varTerm)*dvarTerm;
        else
            dL_dK = dL_dK - ((floor(numData/2)+1)/varTerm)*dvarTerm;
        end
    end
else
    Yscaled = Y;
    for d=1:D
        Yscaled(:,d) = w(d)*Y(:,d);
    end

    dL_dK = -D/2*invK + .5*invK*Yscaled*Yscaled'*invK;
end

dL_dx = zeros(N, q);

for d = 1:q
    Kpart = kernelDiffX(X, theta, d, cacheK);
    dL_dx(:, d) = 2*sum(dL_dK.*Kpart, 2) - diag(dL_dK).*diag(Kpart);
end

% Dynamics Part

if (priorModel(3) == -1)
    dLp_dx = -X;
else
    [Kp, invKp, cacheKp] = computePriorKernel(Xin, thetap, priorModel(3));
    if MARGINAL_DW == 1
        Nq = Np*q;
        dLp_dK = -q/2*invKp;
        varTerm = 0;
        dvarTerm = 0;
        for d = 1:q
            varTerm = varTerm + 0.5*Xout(:, d)'*invKp*Xout(:, d);
            dvarTerm = dvarTerm - 0.5*Xout(:,d)*Xout(:,d)';
        end
        varTerm = varTerm + 0.5/W_VARIANCE;
        dvarTerm = invKp*dvarTerm*invKp;

        if mod(Nq,2) == 0
            dLp_dK = dLp_dK - ((Nq/2)/varTerm + 0.5/varTerm)*dvarTerm;
        else
            dLp_dK = dLp_dK - ((floor(Nq/2)+1)/varTerm)*dvarTerm;
        end
    else
        dLp_dK = -q/2*invKp + 0.5*invKp*Xout*Xout'*invKp;
    end
    
    dLp_dxin = zeros(Np, qp);
    
    for d = 1:qp
        if (priorModel(3) == 0)
            Kpart = kernel2DiffX(Xin, thetap, d);
        elseif (priorModel(3) == 1)
            Kpart = kernelDiffX(Xin, thetap, d, cacheKp);
        elseif (priorModel(3) == 2)
            Kpart = lin_kernelDiffX(Xin, thetap, d);
        elseif (priorModel(3) == 3)
            Kpart = lin_kernel2DiffX(Xin, thetap(1:2), d) + ...
                kernel2DiffX(Xin, thetap(3:5), d);
        elseif (priorModel(3) == 4)
            Kpart = lin_kernel2DiffX(Xin, thetap, d);
        elseif (priorModel(3) == 5)
            Kpart = lin_kernelDiffX(Xin, thetap(1), d) + ...
                kernelDiffX(Xin, thetap(2:3), d, cacheKp);
        end
        dLp_dxin(:, d) = ...
            2*sum(dLp_dK.*Kpart, 2) - diag(dLp_dK).*diag(Kpart);
    end

    if MARGINAL_DW == 1
        dvarTerm = -invKp*Xout;

        if mod(Nq,2) == 0
            dLp_dxout = ((Nq/2)/varTerm + 0.5/varTerm)*dvarTerm;
        else
            dLp_dxout = ((floor(Nq/2)+1)/varTerm)*dvarTerm;
        end
        dLp_dx = priorDiffX(dLp_dxin, dLp_dxout, N, q, seg, priorModel);
    else
        dLp_dx = priorDiffX(dLp_dxin, -invKp*Xout, N, q, seg, priorModel);
    end

    if (priorModel(1) == 0 || priorModel(1) == 1)
        dLp_dx(seg+1,:) = dLp_dx(seg+1,:) - (X(seg+1,:) - X(seg,:));
        dLp_dx(seg,:) = dLp_dx(seg,:) - X(seg,:) + (X(seg+1,:) - X(seg,:));
    elseif (priorModel(1) == 2 || priorModel(1) == 3)
        dLp_dx(seg,:) = dLp_dx(seg,:) - X(seg,:);
    end
end

dL_dx = dL_dx + dLp_dx;
gX= dL_dx(:)';
g = -[gX(:)'];