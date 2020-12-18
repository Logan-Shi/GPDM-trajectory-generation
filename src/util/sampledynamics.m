function X_samples  = sampledynamics(inputModel, oldX_samples, output, hmcoptions)

load(inputModel);
[Xin Xout] = priorIO(X, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
q = size(X,2);
if (modelType(1) < 2)
    M = size(oldX_samples,1) - 2;
    newX = oldX_samples(3:end,:);
    initX = oldX_samples(1:2,:);
else
    M = size(oldX_samples,1) - 1;
    newX = oldX_samples(2:end,:);
    initX = oldX_samples(1,:);
end

[samples] = hmc_wrapper('gpdmdynlikelihoodX', newX(:)', hmcoptions, 'gpdmdyngradientX', ...
    X, thetad, Kd, invKd, segments,  modelType, initX);

for iters = 1:size(samples,1)
    newX = reshape(samples(iters,1:size(samples,2)), size(samples,2)/q, q);
    X_samples{iters} = oldX_samples;
    if (modelType(1) < 2)
        X_samples{iters}(1:2,:) = initX;
        X_samples{iters}(3:end,:) = newX;
    else
        X_samples{iters}(1,:) = initX;
        X_samples{iters}(2:end,:) = newX;
    end
end

save(output, 'X_samples');