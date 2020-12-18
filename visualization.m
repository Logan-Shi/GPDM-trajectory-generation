load example_model
[K invK] = computeKernel(X, theta);
[Xin Xout] = priorIO(X, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
% starts at ened of training sequence;
for plotType = 0:2
    figure()
    gpdmvisualise(X, Y, segments, invK, theta, invKd, thetad,  modelType, plotType)
end