function L = gpdmincompletelikelihood2(params, Y, w, segments, priorModel, thetaSamples)

R = size(thetaSamples,1);

L = 0;
for r = 1:R
    thetaR = exp(thetaSamples(r,1:3)); 
    thetapR = exp(thetaSamples(r,4:end)); 
    
    if (r > 1 && max(thetaSamples(r-1,:) - thetaSamples(r,:)) == 0)
        curL = lastL;
    else
        curL = gpdmposterior(params, Y, w, segments, priorModel, [thetaR thetapR]); 
    end

    L = L + curL;
    lastL = curL;
end

L = L/R;