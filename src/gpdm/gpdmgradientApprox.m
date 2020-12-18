function g = gpdmgradientApprox(test_x, test_y, X, Y, w, A, invK, ...
    theta, Xin, Xout, Ad, invKd, thetad, modelType, missing)


q = size(X,2);
D = size(Y,2);
M = length(test_x)/q;

testX = reshape(test_x, M, q);
testY = reshape(test_y, M, D);

dL_dx = zeros(M, q);

if (modelType(1) < 2)
	order = 2;
else
	order = 1;
end
for n=1:M
    x = testX(n, :);
    
    if (~ismember(n, missing))
        y = testY(n, :);
        
        gn = zeros(1,q);
        
        kbold = kernel(x, X, theta)';
        f = A*kbold;
        sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 
        yHat_i = w.*(y' - f);
        
        if sigma2 < 1e-6;
            sigma2 = 1e-6;
        end
        
        prePart = ...
		theta(1)/sigma2*(w'.*yHat_i'*Y'+(D-yHat_i'*yHat_i/sigma2)*kbold')*invK;
        
		for k = 1:q
            gn(k) = prePart*((x(k) - X(:, k)).*kbold);
        end
    else
		gn = zeros(1,q);
	end

	% dynamic part
    
    if (n < M-(order-1))
        next_x = testX(n+1,:);
        
        if (modelType(1) == 0)
            next2_x = testX(n+2,:);
            input = [next_x x];
            output = next2_x;
        elseif (modelType(1) == 1)
            next2_x = testX(n+2,:);
            input = [next_x next_x - x];
            output = next2_x;
        elseif (modelType(1) == 2)
            input = [x];
            output = next_x;
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
        
        if (modelType(2) == 0)
            xHat_i = (output' - f); 
        elseif (modelType(2) == 1)
            if (modelType(1) < 2)
                xHat_i = (output' - (f + next_x'));
                dL_dx(n+1,:) = dL_dx(n+1,:) - xHat_i'/sigma2;
            else
                xHat_i = (output' - (f + x'));
                dL_dx(n,:) = dL_dx(n,:) - xHat_i'/sigma2;
            end
        end
        
        if sigma2 < 1e-6;
            sigma2 = 1e-6;
        end
        
        if (modelType(1) < 2)
            % specific to 2nd order models
            if (modelType(3) == 0)
                prePart = (xHat_i'*Xout'+(q-xHat_i'*xHat_i/sigma2)*kbold')*invKd;
                for k = 1:2*q
                    if (k < q+1)
                        dL_dx(n+1,k) = dL_dx(n+1,k) +...
                            thetad(1)/sigma2*prePart*((next_x(k) - Xin(:, k)).*kbold);  
                    else
                        if (modelType(1) == 0)
                            dL_dx(n,k-q) = dL_dx(n,k-q) + ...
                                thetad(2)/sigma2*prePart*((x(k-q) - Xin(:, k)).*kbold);
                        elseif (modelType(1) == 1)
                            dL_dx(n,k-q) = dL_dx(n,k-q) - ...
                                thetad(2)/sigma2*prePart*((next_x(k-q) - x(k-q) - Xin(:, k)).*kbold);
                            dL_dx(n+1,k-q) = dL_dx(n+1,k-q) + ...
                                thetad(2)/sigma2*prePart*((next_x(k-q) - x(k-q) - Xin(:, k)).*kbold);
                        end
                    end
                end
            elseif (modelType(3) == 1)
                fprintf('ERROR: bad modelType');
            elseif (modelType(3) == 2)
                
            end
        elseif (modelType(1) == 2)
            % first order
            if (modelType(3) == 0)
                fprintf('ERROR: bad modelType');
            elseif (modelType(3) == 1)
                prePart = thetad(1)/sigma2*(xHat_i'*Xout'+(q-xHat_i'*xHat_i/sigma2)*kbold')*invKd;
                for k = 1:q
                    dL_dx(n,k) = dL_dx(n,k) + prePart*((x(k) - Xin(:, k)).*kbold);
                end
            elseif (modelType(3) == 2)
                % this doesn't work 
                prePart = thetad(1)/sigma2*(xHat_i'*Xout'+(q-xHat_i'*xHat_i/sigma2)*kbold')*invKd;
                for k = 1:q
                    dL_dx(n,k) = dL_dx(n,k) + prePart*Xin(:, k);
                end
                dL_dx(n,:) = dL_dx(n,:) + thetad(1)*((q-xHat_i'*xHat_i/sigma2)/(sigma2))*x;
            end
        end
    end
   
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
                
        if (modelType(2) == 0)
            xHat_i = (x' - f); 
        elseif (modelType(2) == 1)
            xHat_i = (x' - (f + last_x'));
        end
        
        gn = gn + ((xHat_i)/sigma2)';
    end
	
	if (n == 2 && order == 2)
        last_x = testX(n-1,:);
        gn = gn + x - last_x;
	end
    
	if (n == 1)
		if (order == 2)
			next_x = testX(n+1,:);
        	gn = gn + x - (next_x - x);
		elseif (order == 1)
			gn = gn + x;
		end
    end

    dL_dx(n,:) = dL_dx(n,:) + gn;
end

g = [dL_dx(:)'];




