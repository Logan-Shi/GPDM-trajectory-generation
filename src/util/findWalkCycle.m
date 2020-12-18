function [Out] = findWalkCycle(Y)

[X, I] = sort(Y(:,59) + Y(:,52));
Out = I(1);
for n=2:length(I)
    if (min(abs(I(n) - Out)) > 10)
        Out = [Out I(n)];
    end
end

