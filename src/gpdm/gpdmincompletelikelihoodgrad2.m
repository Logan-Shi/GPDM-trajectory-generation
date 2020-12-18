function g = gpdmincompletelikelihoodgrad2(params, Y, w, segments, priorModel, thetaSamples)

R = size(thetaSamples,1);

g = zeros(size(params));
for r = 1:R
    thetaR = exp(thetaSamples(r,1:3)); 
    thetapR = exp(thetaSamples(r,4:end)); 
    
    if (r > 1 && max(thetaSamples(r-1,:) - thetaSamples(r,:)) == 0)
        curG = lastG;
    else
        curG = gpdmposteriorgrad(params, Y, w, segments, priorModel, [thetaR thetapR]);
    end

    g = g + curG;
    lastG = curG;
end

g = g/R;
