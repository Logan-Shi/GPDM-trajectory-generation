function X_opt  = optimizedynamics(inputModel, oldX_opt, output, options, extIters)

load(inputModel);
[Xin Xout] = priorIO(X, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
q = size(X,2);
if (modelType(1) < 2)
    M = size(oldX_opt,1) - 2;
    newX = oldX_opt(3:end,:);
    initX = oldX_opt(1:2,:);
else
    M = size(oldX_opt,1) - 1;
    newX = oldX_opt(2:end,:);
    initX = oldX_opt(1,:);
end

for iters = 1:extIters
    fprintf(2,'Iteration %d\n',iters);

    params = [newX(:)'];
    params = scg('gpdmdynlikelihoodX', params, options, 'gpdmdyngradientX',...
        X, thetad, Kd, invKd, segments, modelType, initX);
    gpdmdynlikelihoodX(params, X, thetad, Kd, invKd, segments, modelType, initX)
    newX = reshape(params(1:M*q), M, q);
end

X_opt = oldX_opt;
if (modelType(1) < 2)
    X_opt(1:2,:) = initX;
    X_opt(3:end,:) = newX;
else
    X_opt(1,:) = initX;
    X_opt(2:end,:) = newX;
end

save(output, 'X_opt');



