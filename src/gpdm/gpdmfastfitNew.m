function [Xt] = gpdmfastfitNew(Xt, Yt, X, Y, w, theta, invK, ...
    Xin, Xout, invKd, thetad, opt, extItr, modelType, missing)

% This function only works for 2 new elements.  It might be fixed to work 
% in general later.
	
if (~exist('missing', 'var'))
    missing = [];
end


N = size(Y,1);
D = size(Y,2);
q = size(X,2);
M = size(Yt,1);
test_y = [Yt(:)'];
A = Y'*invK;
Ad = Xout'*invKd;

for iters = 1:extItr
    %fprintf(2,'Iteration %d\n',iters);
    
    test_x = [Xt(:)'];
    %gpdmtestlikelihood(test_x, test_y, X, Y, w, A, ...
     %   invK, theta, Xin, Xout, Ad, invKd, thetad, modelType, missing)
    test_x = scg('gpdmlikelihoodApprox', test_x, opt, ...
        'gpdmgradientApprox', test_y, X, Y, w, A, ...
        invK, theta, Xin, Xout, Ad, invKd, thetad, modelType, missing);

    Xt = reshape(test_x(1:M*q), M, q);
end
				
