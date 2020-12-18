load example_model
[K invK] = computeKernel(X, theta);
[Xin Xout] = priorIO(X, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
simSteps = 256;
% starts at ened of training sequence;
simStart = [X(segments(1)+1,:) X(end,:)]; %  inputs 2 points in case using 2nd order model
[X_pred XRand_pred] = simulatedynamics(X, segments, thetad, invKd, simSteps, simStart, modelType);

% uncomment if want to generate new samples

hmcopt = foptions;      % Default options vector.
hmcopt(1) = 1;			% Switch on diagnostics.
hmcopt(7) = 100;	    	% Number of steps in leapfrog trajectory.
hmcopt(9) = 0; 
hmcopt(14) = 60;		% Number of Monte Carlo samples returned. 
hmcopt(15) = 20;		% Number of samples omitted at start of chain.
hmcopt(18) = 0.01;  	% leapfrog step size
X_samples = sampledynamics('example_model', X_pred, 'samples', hmcopt);

% load samples
% clf;
% hold on;

for n=1:size(X_samples,2)
plotseries(X_samples{n}, [1], 'g');
end
plotseries(X, segments, 'b');
plotseries(X_pred, [1], 'r');
