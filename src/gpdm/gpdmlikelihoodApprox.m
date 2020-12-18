function L = gpdmlikelihoodApprox(test_x, test_y, X, Y, w, A, ...
    invK, theta, Xin, Xout, Ad, invKd, thetad, modelType, missing)
% A = Y'*invK
% Ad = Xout'*invKd

q = size(X,2);
D = size(Y,2);
M = length(test_x)/q;
testX = reshape(test_x, M, q);
testY = reshape(test_y, M, D);
L = 0;
if (modelType(1) < 2)
	order = 2;
else
	order = 1;
end

for n=1:M
    x = testX(n,:);
    if (~ismember(n, missing))
        y = testY(n,:);
        
        kbold = kernel(x, X, theta)';
        f = A*kbold;
        sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 
        
        if sigma2 < 1e-6;
            sigma2 = 1e-6;
        end
        
        yHat_i = w.*(y' - f);
        
        L = L -D/2*log(sigma2) -D/2*log(2*pi) - yHat_i'*yHat_i/(2*sigma2);
    end
    
    % dynamic part
	if (n > order)
        last_x = testX(n-1,:);
		if (modelType(1) == 0)
        	last2_x = testX(n-2,:);
			input = [last_x last2_x];
		elseif (modelType(1) == 1)
        	last2_x = testX(n-2,:);
			input = [last_x last_x - last2_x];
		elseif (modelType(1) == 2)
			input = [last_x];
		end
        
        if (modelType(3) == 0)
            kbold = kernel2(input, Xin, thetad)';
			sigma2 = thetad(3) +1/thetad(end) - kbold'*invKd*kbold;
        elseif (modelType(3) == 1)
            kbold = kernel(input, Xin, thetad)';
			sigma2 = thetad(2) +1/thetad(end) - kbold'*invKd*kbold;
        elseif (modelType(3) == 2)
            kbold = lin_kernel(input, Xin, thetad)';
            sigma2 = lin_kernel(input, input, thetad) +1/thetad(end) - kbold'*invKd*kbold;
        end

        f = Ad*kbold;
        if sigma2 < 1e-6;
            sigma2 = 1e-6;
        end
        
        if (modelType(2) == 0)
            xHat_i = (x' - f);  
        elseif (modelType(2) == 1)
            xHat_i = (x' - (f + last_x'));
        end
        
        L = L -q/2*log(sigma2) -q/2*log(2*pi) - xHat_i'*xHat_i/(2*sigma2);
    end
	
	if (n == 2 && order == 2)
        last_x = testX(n-1,:);
        L = L - (x - last_x)*(x - last_x)'/2;
	end
    
	if (n == 1)
        L = L - x*x'/2;
    end
end

L = -L;

