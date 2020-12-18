function g = gpdmincompletelikelihoodgrad(params, Y, w, segments, priorModel, Xsamples, q)

R = size(Xsamples,1);

g = zeros(size(params));
for r = 1:R
    Xr = reshape(Xsamples(r,1:size(Xsamples,2)), size(Xsamples,2)/q, q);

    if (r > 1 && max(Xsamples(r-1,:) - Xsamples(r,:)) == 0)
        curG = lastG;
    else
        curG = gpdmgradient(params, Y, w, segments, priorModel, [], Xr);
    end

    g = g + curG;
    lastG = curG;
end

g = g/R;


