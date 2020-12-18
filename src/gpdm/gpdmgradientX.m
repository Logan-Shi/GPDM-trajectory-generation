function g = gpdmgradientX(params, Y, w, segments, priorModel, ...
missing, fixedX, fixedTheta)

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

if (~exist('fixedX','var'))
    q = (length(params)-(ndp+3))/N;
    X = reshape(params(1:N*q), N, q);
    M = N;
    seg = segments;
    [Xin Xout] = priorIO(X, seg, priorModel);
    Np = size(Xin, 1);
    Mp = Np;
else
    q = size(fixedX, 2);
    M = N - size(fixedX, 1);
    X = [fixedX; reshape(params(1:M*q), M, q)];
    if (exist('fixedTheta', 'var'))
        seg = [segments N-M+1];
    else
        seg = segments;
    end
    [Xin Xout] = priorIO(X, seg, priorModel);
    Np = size(Xin, 1);
	if (priorModel(1) == 0 || priorModel(1) == 1)
    	Mp = M - 2;
	else
		Mp = M - 1;
	end
end

if (~exist('missing', 'var'))
    nmissing = 1:N;
else
    if (exist('fixedX', 'var') && (exist('fixedTheta', 'var')))
        nmissing = setdiff(1:N, missing+N-M);
    else
        nmissing = setdiff(1:N, missing);
    end
end

qp = size(Xin, 2);

if (~exist('fixedTheta', 'var'))
    lntheta = params(end-(ndp+2):end-ndp);
    lnthetap = params(end-(ndp-1):end);
    theta = exp(lntheta);
    thetap = exp(lnthetap);
    theta = thetaConstrain(theta);
    thetap = thetaConstrain(thetap);
else
    theta = fixedTheta(end-(ndp+2):end-ndp);
    thetap = fixedTheta(end-(ndp-1):end);
end

[K, invK] = computeKernel(X(nmissing,:), theta);


Yscaled = Y(nmissing,:);
for d=1:D
    Yscaled(:,d) = w(d)*Y(nmissing,d); 
end

dL_dK = -D/2*invK + .5*invK*Yscaled*Yscaled'*invK;

if (~exist('fixedTheta', 'var'))
    dk = zeros(1, 3);
    [dK{1}, dK{2}] = kernelDiffParams(X(nmissing,:), theta);
    for i = 1:2
        dk(i) = sum(sum(dL_dK.*dK{i}));
    end
    dk(3) = -sum(diag(dL_dK)/theta(end).^2);
    
    gParam = dk.*theta-1;
end

dL_dx = zeros(N, q);

for d = 1:q
    Kpart = kernelDiffX(X(nmissing,:), theta, d);
    dL_dx(nmissing, d) = 2*sum(dL_dK.*Kpart, 2) - diag(dL_dK).*diag(Kpart);
end

% Dynamics Part

    gParamp = zeros(1,ndp);
    
    
    if (priorModel(3) == -1)
        dLp_dx = -X;
    else
        [Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));
        dLp_dK = -q/2*invKp + .5*invKp*Xout*Xout'*invKp;
        if (~exist('fixedTheta', 'var'))
            
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
        end
        
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
    end
    
dL_dx = dL_dx + dLp_dx;

gParam = [0 0 0];
gParamp = zeros(1, ndp);
%dL_dx(missing,:) = zeros(length(missing), q);
if (~exist('fixedTheta', 'var') && ~exist('fixedX', 'var') )
    gX= dL_dx(:)';
    g = -[gX(:)' gParam gParamp];
elseif (exist('fixedX', 'var') && ~exist('fixedTheta','var'))
    g = -[gParam gParamp];
else
    gX = dL_dx(N-M+1:end,:);
    g = -[gX(:)'];
end

