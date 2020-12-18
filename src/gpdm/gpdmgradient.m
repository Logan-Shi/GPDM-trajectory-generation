function g = gpdmgradient(params, Y, w, segments, priorModel, ...
missing, fixedX, fixedTheta, Kn, invKn, Kdn, invKdn)

global USE_GAMMA_PRIOR
global GAMMA_ALPHA
global LAWRENCE_PARAMS
global USE_LAWRENCE
global FIX_HP
global MARGINAL_W
global MARGINAL_DW
global REMOVE_REDUNDANT_SCALE
global W_VARIANCE
global M_CONST
global BALANCE
global USE_OLD_MISSING_DATA

% global LOGTHETA

% Could encourage log(theta) to be close to log(thetad).
tieHps = 0;
if (tieHps == 1)
power = 1;
mult = 0;
end

recConst = M_CONST;
lambda = BALANCE; 

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
    missing = []; 
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
    
%     if ~LOGTHETA
        theta = exp(lntheta);
        thetap = exp(lnthetap); 
%     else
%         theta = [lntheta(1:end-1) exp(lntheta(end))];
%         thetap = [lnthetap(1:end-1) exp(lnthetap(end))];
%     end
else
    theta = fixedTheta(end-(ndp+2):end-ndp);
    thetap = fixedTheta(end-(ndp-1):end);
end
% theta = thetaConstrain(theta);
% thetap = thetaConstrain(thetap);

reuse = 0; 
if (exist('invKdn', 'var'))
    reuse = 1; 
end

if (reuse == 1)
    [K, invK] = computeKernel(X(nmissing,:), theta, Kn, invKn);
else
    [K, invK] = computeKernel(X(nmissing,:), theta);
end

if (USE_OLD_MISSING_DATA == 1) 
else
if (~isempty(missing))
    kbold = kernel(X(missing+N-M,:), X(nmissing,:), theta)';
    A = Y(nmissing,:)'*invK*kbold;
    Y(missing+N-M,:) = A';
    nmissing = 1:N;
    if (reuse == 1)
        [K, invK] = computeKernel(X(nmissing,:), theta, Kn, invKn);
    else
        [K, invK] = computeKernel(X(nmissing,:), theta);
    end
end
end
cacheK = K - eye(size(X(nmissing,:), 1))*1/theta(end);

numData = length(nmissing); 
if MARGINAL_W == 1
    dL_dK = -D/2*invK;
    for d = 1:D
        varTerm = 0.5*Y(nmissing, d)'*invK*Y(nmissing, d)+0.5/W_VARIANCE;
        dvarTerm = -0.5*invK*Y(nmissing, d)*Y(nmissing, d)'*invK; 
        if mod(numData,2) == 0
            dL_dK = dL_dK - ((numData/2)/varTerm + 0.5/varTerm)*dvarTerm; 
        else
            dL_dK = dL_dK - ((floor(numData/2)+1)/varTerm)*dvarTerm; 
        end
    end 
else
    Yscaled = Y(nmissing,:);
    for d=1:D
        Yscaled(:,d) = w(d)*Y(nmissing,d);
    end
    dL_dK = -D/2*invK + .5*invK*Yscaled*Yscaled'*invK;
end

gParam = zeros(size(theta)); 
gParamp = zeros(size(thetap)); 

if (~exist('fixedTheta', 'var'))
    dk = zeros(1, 3);
    [dK{1}, dK{2}] = kernelDiffParams(X(nmissing,:), X(nmissing,:), theta, cacheK);

    for i = 1:2
        dk(i) = sum(sum(dL_dK.*dK{i}));
    end
    dk(3) = -sum(diag(dL_dK)/theta(end).^2);
    
    %
    if (tieHps == 1)
        gParam = dk.*theta+ ... 
            2*power*exp(mult)*prod(theta)/(prod(thetap)-exp(mult)*prod(theta))-recConst;
    else
%         if ~LOGTHETA
            %gParam = dk.*theta; 
            gParam = dk.*theta - recConst; 
%         else
%             gParam(1:end-1) = (dk(1:end-1).*theta(1:end-1))./theta(1:end-1) - ...
%                 recConst./theta(1:end-1);
%             gParam(end) = dk(end)*theta(end) - recConst;
%         end
    end
end

dL_dx = zeros(N, q);

for d = 1:q
    Kpart = kernelDiffX(X(nmissing,:), theta, d, cacheK);
    dL_dx(nmissing, d) = 2*sum(dL_dK.*Kpart, 2) - diag(dL_dK).*diag(Kpart);
end

% Dynamics Part

gParamp = zeros(1,ndp);

if (priorModel(3) == -1)
    dLp_dx = -X;
else
    if (reuse == 1) 
        [Kp, invKp, cacheKp] = computePriorKernel(Xin, thetap, priorModel(3), Kdn, invKdn);
    else
        [Kp, invKp, cacheKp] = computePriorKernel(Xin, thetap, priorModel(3));
    end

    
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
        elseif (USE_GAMMA_PRIOR == 1)
            [scale snr_2 invl_2 dscale dsnr dinvl] = alpha_info(thetap, priorModel);
            %[G, dG_dscale] = gammaPrior(GAMMA_ALPHA(1), LAWRENCE_PARAMS(1), scale); 
            [G, dG_dsnr_2] = gammaPrior(GAMMA_ALPHA(2), LAWRENCE_PARAMS(2), snr_2); 
            [G, dG_dinvl_2] = gammaPrior(GAMMA_ALPHA(3), LAWRENCE_PARAMS(3), invl_2); 
            %gParamp = dk.*thetap + (dG_dscale*dscale + dG_dsnr_2*dsnr + ... 
             %   dG_dinvl_2*dinvl).*thetap;  

             gParamp = dk.*thetap + (dG_dsnr_2*dsnr + ... 
                dG_dinvl_2*dinvl).*thetap;  
%              [G, dG_dt13] = gammaPrior(GAMMA_ALPHA(1), LAWRENCE_PARAMS(1), thetap(1) + thetap(3)); 
%              [G, dG_dt2] = gammaPrior(GAMMA_ALPHA(2), LAWRENCE_PARAMS(2), thetap(2)); 
%              [G, dG_dt4] = gammaPrior(GAMMA_ALPHA(3), LAWRENCE_PARAMS(3), thetap(4)); 
% 
% 
%              gParamp = dk.*thetap + (dG_dt13*[1 0 1 0] + dG_dt2*[0 1 0 0] + ... 
%                  dG_dt4*[0 0 0 1]).*thetap; 
%             phi = 1./(thetap.*thetap);
%             phi(end) = 1/phi(end);
%             omega(2) = thetap(3)*thetap(3)/100^4;
%             %omega(2) - phi(3)/100^4
%             priorTmp = (0.5*alpha./phi - 1./phi - 0.5*alpha./omega);
%             priorTmp(end) = priorTmp(end) + (0.5*100^4*phi(2)*alpha(2)/(phi(3)*phi(3)) - 0.5*alpha(2)/phi(3));
%             priorTmp(1:end-1) = priorTmp(1:end-1).*(-2.*(thetap(1:end-1).^(-3))).*thetap(1:end-1);
%             priorTmp(end) = priorTmp(end).*(2*thetap(end)).*thetap(end);
%             gParamp = dk.*thetap + priorTmp;
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

%gParam = [0 0 0]; 

if (~exist('fixedTheta', 'var') && ~exist('fixedX', 'var') )
    gX= dL_dx(:)';
    g = -[gX(:)' gParam gParamp];
elseif (exist('fixedX', 'var') && ~exist('fixedTheta','var'))
    g = -[gParam gParamp];
else
    gX = dL_dx(N-M+1:end,:);
    g = -[gX(:)'];
end

