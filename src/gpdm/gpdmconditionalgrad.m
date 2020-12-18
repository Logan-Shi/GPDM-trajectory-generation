function g = gpdmconditionalgrad(params, Y, w, segments, priorModel, ...
    missing, Xn, fixedTheta, Kn, invKn, Kdn, invKdn)

% TODO: still have to work this out.

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
cacheB = B - eye(size(Xm(nmissing,:), 1))*1/theta(end);

Km = B - A'*invKn*A; 
invKm = pdinv(Km); 

    [Apart Bpart] = condKernelDiffX(Xm, Xn, theta, 1, A, cacheB);
    dA_dx = Apart; 
    dA_dx(:,2:end) = 0;
    dB_dx = Bpart; 
    dB_dx(2:end, 2:end) = 0; 
     dC_dx = dB_dx - (dA_dx'*invKn*A + A'*invKn*dA_dx); 

numData = length(nmissing); 

Yscaled = Ym(nmissing,:);
dL_dx = zeros(M, q);
dL_dx(1,1) = -D/2*trace(invKm*dC_dx); 
for d=1:D
     Yscaled(:,d) = w(d)*(Ym(nmissing, d) - A'*invKn*Yn(:,d)); 
     dmu_dx = dA_dx'*invKn*w(d)*Yn(:,d); 
     dL_dx(1,1) = dL_dx(1,1) + .5*Yscaled(:,d)'*invKm*dC_dx*invKm*Yscaled(:,d) + Yscaled(:,d)'*invKm*dmu_dx; 
end

gX = dL_dx;
g = -[gX(:)'];

return; 
% Dynamics Part

gParamp = zeros(1,ndp);

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
        dLp_dK = lambda*dLp_dK; 
    end

    if (~exist('fixedTheta', 'var'))

        dk = zeros(1, ndp);
        if (priorModel(3) == 0)
            [dK{1}, dK{2}, dK{3}] = kernel2DiffParams(Xin, thetap);
        elseif (priorModel(3) == 1)
            [dK{1}, dK{2}] = kernelDiffParams(Xin, Xin, thetap, cacheKp);
        elseif (priorModel(3) == 2)
            [dK{1}] = lin_kernelDiffParams(Xin, thetap);
        elseif (priorModel(3) == 3)
            [dK{1} dK{2}] = lin_kernel2DiffParams(Xin, thetap(1:2));
            [dK{3}, dK{4}, dK{5}] = kernel2DiffParams(Xin, thetap(3:5));
        elseif (priorModel(3) == 4)
            [dK{1} dK{2}] = lin_kernel2DiffParams(Xin, thetap);
        elseif (priorModel(3) == 5)
            [dK{1}] = lin_kernelDiffParams(Xin, thetap(1));
            [dK{2}, dK{3}] = kernelDiffParams(Xin, Xin, thetap(2:3), cacheKp);
        end

        for i = 1:(ndp-1)
            dk(i) = sum(sum(dLp_dK.*dK{i}));
        end
        dk(ndp) = -sum(diag(dLp_dK)/thetap(end).^2);


        %
        if (tieHps == 1)
            gParamp = dk.*thetap...
                -2*power*prod(thetap)/(prod(thetap)-exp(mult)*prod(theta));
        elseif (gammaPrior == 1)
            phi = 1./(thetap.*thetap);
            phi(end) = 1/phi(end);
            omega(2) = thetap(3)*thetap(3)/100^4;
            %omega(2) - phi(3)/100^4
            priorTmp = (0.5*alpha./phi - 1./phi - 0.5*alpha./omega);
            priorTmp(end) = priorTmp(end) + (0.5*100^4*phi(2)*alpha(2)/(phi(3)*phi(3)) - 0.5*alpha(2)/phi(3));
            priorTmp(1:end-1) = priorTmp(1:end-1).*(-2.*(thetap(1:end-1).^(-3))).*thetap(1:end-1);
            priorTmp(end) = priorTmp(end).*(2*thetap(end)).*thetap(end);
            gParamp = dk.*thetap + priorTmp;
        else
            %                 if ~LOGTHETA
            %gParamp = dk.*thetap;
            gParamp = dk.*thetap - 1;
            %                     gParamp = dk.*thetap;
            %                     gParamp(1) = gParamp(1) - (thetap(1) - 0.2)*thetap(1);
            %                     gParamp(2) = gParamp(2) - (thetap(2) - 10000/thetap(3))*thetap(2);
            %                     gParamp(3) = gParamp(3) - (thetap(2) - 10000/thetap(3))*10000/thetap(3) - 1;
            %                 else
            %                     gParamp(1:end-1) = (dk(1:end-1).*thetap(1:end-1))./thetap(1:end-1) - ...
            %                         1./thetap(1:end-1);
            %                     gParamp(end) = dk(end)*thetap(end) - 1;
            %                 end
        end
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
        dLp_dx = priorDiffX(dLp_dxin, -lambda*invKp*Xout, N, q, seg, priorModel);
    end

    if (priorModel(1) == 0 || priorModel(1) == 1)
        dLp_dx(seg+1,:) = dLp_dx(seg+1,:) - lambda*(X(seg+1,:) - X(seg,:));
        dLp_dx(seg,:) = dLp_dx(seg,:) - lambda*X(seg,:) + lambda*(X(seg+1,:) - X(seg,:));
    elseif (priorModel(1) == 2 || priorModel(1) == 3)
        dLp_dx(seg,:) = dLp_dx(seg,:) - lambda*X(seg,:);
    end
end

dL_dx = dL_dx + dLp_dx;
%dL_dx(seg,:) = 0;

%norm(dL_dx)
%norm(dLp_dx)
%gParam = [0 0 0];
%gParamp = [0 0 0];
%dL_dx(missing,:) = zeros(length(missing), q);
%gParam(end) = 0;
%gParamp(1) = 0;
% 
% [theta, gParam] = thetaConstrain(theta, gParam);
% [thetap, gParamp] = thetaConstrain(thetap, gParamp);

if USE_LAWRENCE == 1
    gParamp = zeros(1,ndp);
end

if FIX_HP == 1
    gParam = [0 0 0];
    gParamp = zeros(size(gParamp));
end

if REMOVE_REDUNDANT_SCALE == 1
    gParam(2) = 0;
end

if MARGINAL_DW == 1
    if (priorModel(3) == 1)
        gParamp(2) = 0; 
    end
end

if (~exist('fixedTheta', 'var') && ~exist('fixedX', 'var') )
    gX= dL_dx(:)';
    g = -[gX(:)' gParam gParamp];
elseif (exist('fixedX', 'var') && ~exist('fixedTheta','var'))
    g = -[gParam gParamp];
else
    gX = dL_dx(N-M+1:end,:);
    g = -[gX(:)'];
end

