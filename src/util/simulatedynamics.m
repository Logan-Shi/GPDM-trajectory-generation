function [X_pred, XRand_pred] = simulatedynamics(X, segments, thetad, invKd, ...,
    simSteps, simStart, modelType, updateKernel)
if (~exist('updateKernel', 'var'))
    updateKernel = 0;
end
[Xin Xout] = priorIO(X, segments, modelType);

q = size(X,2);

X_pred = zeros(simSteps, q);
XRand_pred = zeros(simSteps, q);

if (modelType(1) < 2)
    order = 2;
    X_pred(2,:) = simStart(:,1:q);
    X_pred(1,:) = simStart(:,end-q+1:end);
    XRand_pred(2,:) = simStart(:,1:q);
    XRand_pred(1,:) = simStart(:,end-q+1:end);
else
    order = 1;
    X_pred(1,:) = simStart(:,end-q+1:end);
    XRand_pred(1,:) = simStart(:,end-q+1:end);
end

curX = X;
curRandX = X;
for n = (order+1):simSteps
    if (updateKernel == 1)
        curX = [curX; X_pred(n-1,:)];
        [Xin Xout] = priorIO(curX, segments, modelType);
        [Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
        curRandX = [curRandX; XRand_pred(n-1,:)];
        [RXin RXout] = priorIO(curRandX, segments, modelType);
        [RKd invRKd] = computePriorKernel(RXin, thetad, modelType(3));
    else
        RXin = Xin;
        RXout = Xout;
        invRKd = invKd;
    end
    if (order == 1)
        [X_pred(n,:) var_pred] = ...
            priorManifoldOutputs([X_pred(n-1,:) zeros(1,q)] , Xin, Xout, thetad, invKd, modelType);
        [XRand_pred(n,:) var_pred] = ...
            priorManifoldOutputs([XRand_pred(n-1,:) zeros(1,q)] , RXin, RXout, thetad, invRKd, modelType);
    else
        [X_pred(n,:) var_pred] = ...
            priorManifoldOutputs([X_pred(n-1,:) X_pred(n-2,:)] , Xin, Xout, thetad, invKd, modelType);
        [XRand_pred(n,:) var_pred] = ...
            priorManifoldOutputs([XRand_pred(n-1,:) X_pred(n-2,:)] , RXin, RXout, thetad, invRKd, modelType);
        
    end
    
%     if (updateKernel == 0)
%         XRand_pred(n,:) = X_pred(n,:) + sqrt(var_pred)*randn(1,q);    
%     else
%         XRand_pred(n,:) = XRand_pred(n-1,:) + sqrt(var_pred)*randn(1,q);
%     end
    XRand_pred(n,:) = XRand_pred(n-1,:) + sqrt(var_pred)*randn(1,q);
end