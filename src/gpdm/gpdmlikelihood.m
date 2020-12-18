function L = gpdmlikelihood(params, Y, w, segments, priorModel, ...
    missing, fixedX, fixedTheta, Kn, invKn, Kdn, invKdn)

global INFO
global USE_GAMMA_PRIOR
global GAMMA_ALPHA
global LAWRENCE_PARAMS
global MARGINAL_W
global MARGINAL_DW
global W_VARIANCE
global M_CONST
global BALANCE
global USE_OLD_MISSING_DATA

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

Np = size(Xin, 1);

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


numData = length(nmissing); 
if MARGINAL_W == 1
    L = D*log(2)-D*log(sqrt(W_VARIANCE))-(0.5*numData+0.5)*D*log(2*pi)-D/2*logdet(K);
    for d = 1:D
        varTerm = 0.5*Y(nmissing, d)'*invK*Y(nmissing, d)+0.5/W_VARIANCE;
        if mod(numData,2) == 0
            L = L + (numData/2)*(log(numData-1) - log(2*varTerm)) + 0.5*(log(pi) - log(4*varTerm));
        else
            L = L + floor(numData/2)*(log(numData-1) - log(2*varTerm)) - log(2*varTerm);
        end
    end
    L = L - recConst*sum(log(theta));
else
    L = -D*numData/2*log(2*pi)-D/2*logdet(K);
    for d= 1:D
        L = L -0.5*w(d)*w(d)*Y(nmissing, d)'*invK*Y(nmissing, d);
    end
    if (INFO == 1)
        fprintf('-L reconstruction: %4.2f\n', -L);
    end
    L = L + numData*sum(log(w)); 
    
  
    
    L = L - recConst*sum(log(theta)); 
    if (W_VARIANCE > 0) 
        L = L + D*log(2) - D/2*log(2*pi*W_VARIANCE) - 0.5/W_VARIANCE*sum(w.*w) ;
    end
end


if (priorModel(3) == -1)
    Lp = -N*q/2*log(2*pi) -0.5*sum(sum(X.*X));
else
    if (reuse == 1) 
        [Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3), Kdn, invKdn);
    else
        [Kp, invKp] = computePriorKernel(Xin, thetap, priorModel(3));
    end

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
        Lp = lambda*Lp; 

        %
        if (tieHps == 1)
            Lp = Lp - power*log((prod(thetap)-exp(mult)*prod(theta))^2);
        elseif (USE_GAMMA_PRIOR == 1)
            [scale snr_2 invl_2] = alpha_info(thetap, priorModel);
%             scale
%             gammaPrior(GAMMA_ALPHA(1), LAWRENCE_PARAMS(1), scale)
            Lp = Lp + ... %gammaPrior(GAMMA_ALPHA(1), LAWRENCE_PARAMS(1), scale) + ...
                gammaPrior(GAMMA_ALPHA(2), LAWRENCE_PARAMS(2), snr_2) + ...
                gammaPrior(GAMMA_ALPHA(3), LAWRENCE_PARAMS(3), invl_2);


%             Lp = Lp + gammaPrior(GAMMA_ALPHA(1), LAWRENCE_PARAMS(1), thetap(1) + thetap(3)) + ... 
%                 gammaPrior(GAMMA_ALPHA(2), LAWRENCE_PARAMS(2), thetap(2)) + ...
%                 gammaPrior(GAMMA_ALPHA(3), LAWRENCE_PARAMS(3), thetap(4)); 
%            

%             phi = 1./(thetap.*thetap);
%             phi(end) = 1/phi(end);
%             omega(2) = thetap(3)*thetap(3)/100^4;
%             Lp = Lp + sum(0.5*alpha.*log(0.5*alpha./omega) + ...
%                 (0.5*alpha - 1).*log(phi) - log(gamma(0.5*alpha)) - ...
%                 0.5*(phi.*alpha)./omega);
        else
            Lp = Lp - sum(log(thetap));
            %         Lp = Lp - 0.5*(thetap(1) - 0.2)*(thetap(1) - 0.2) - 0.5*log(2*pi) - ...
            %             0.5*(thetap(2) - 10000/thetap(3))*(thetap(2) - 10000/thetap(3)) - 0.5*log(2*pi) - ...
            %             log(thetap(3));
        end
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

