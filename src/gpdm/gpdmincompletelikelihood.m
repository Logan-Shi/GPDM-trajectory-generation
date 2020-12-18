function L = gpdmincompletelikelihood(params, Y, w, segments, priorModel, Xsamples, q)

R = size(Xsamples,1);

L = 0;
for r = 1:R
    Xr = reshape(Xsamples(r,1:size(Xsamples,2)), size(Xsamples,2)/q, q);

    if (r > 1 && max(Xsamples(r-1,:) - Xsamples(r,:)) == 0)
        curL = lastL;
    else
        curL = gpdmlikelihood(params, Y, w, segments, priorModel, [], Xr);
    end

    L = L + curL;
    lastL = curL;
end

L = L/R;