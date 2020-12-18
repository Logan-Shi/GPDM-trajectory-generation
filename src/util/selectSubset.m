function [I, segments] = selectSubset(old_Y, old_segments, numPoints)

n = ceil(size(old_Y,1)*rand); 

K = old_Y - repmat(old_Y(n,:), size(old_Y,1), 1); 
K = sum(K.*K, 2); 
[K I] = sort(K); 

I = I(1:numPoints); 
I = sort(I); 
segments = [1]; 
for k=2:size(I)
    if (sum(ismember(old_segments, I(k))) > 0)
        segments = [segments k]; 
    elseif (I(k) ~= (I(k-1) + 1))
        segments = [segments k]; 
    end
end

