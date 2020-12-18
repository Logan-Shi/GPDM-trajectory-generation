function [Xin, Xout] = priorIO(X, segments, priorModel) 

q = size(X,2);

if (priorModel(1) == 0)
    Xin = [[zeros(2,q); X(2:end-1,:)] [zeros(2,q); X(1:end-2,:)]];
    Xin(union(segments,segments+1),:) = [];
elseif (priorModel(1) == 1)
    Xin = [[zeros(2,q); X(2:end-1,:)] [zeros(2,q); ...
	X(2:end-1,:) - X(1:end-2,:)]];
    Xin(union(segments,segments+1),:) = [];
elseif (priorModel(1) == 2)
	Xin = [zeros(1,q); X(1:end-1,:)];
	Xin(segments,:) = [];
elseif (priorModel(1) == 3)
    Xin = [zeros(1,q); X(1:end-1,:)];
    Xin(segments,:) = [];
    Xin = [Xin ones(size(Xin,1), 1)];
end

if (priorModel(2) == 0)
    Xout = X;
	if (priorModel(1) == 2 || priorModel(1) == 3)
		Xout(segments,:) = [];
	else
    	Xout(union(segments,segments+1),:) = [];
	end
elseif (priorModel(2) == 1)
    Xout = [zeros(1,q); X(2:end,:) - X(1:end-1,:)];
	if (priorModel(1) == 2 || priorModel(1) == 3)
		Xout(segments,:) = [];
	else
    	Xout(union(segments,segments+1),:) = [];
	end
end
