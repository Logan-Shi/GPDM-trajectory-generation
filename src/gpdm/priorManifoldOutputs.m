function [xo, varxo] = priorManifoldOutputs(xi, Xin, Xout, thetad, invKd, modelType)

q = size(Xout,2);

if (modelType(1) == 0)
    input = xi;
elseif (modelType(1) == 1)
    input = [xi(:,1:q) xi(:,1:q)-xi(:,q+1:end)];
elseif (modelType(1) == 2)
    input = [xi(:,1:q)];
end

if (modelType(3) == 0)
    [output, output_var] = ...
        manifoldOutputs2(input, Xin, Xout, thetad, invKd);
elseif (modelType(3) == 1)
    [output, output_var] = ...
        manifoldOutputs(input, Xin, Xout, thetad, invKd);
elseif (modelType(3) == 2)
    [output, output_var] = ...
        linManifoldOutputs(input, Xin, Xout, thetad, invKd);
elseif (modelType(3) == 3)
    N = size(Xout, 1);
    M = size(input, 1);
    
    alpha = zeros(N, q);
    
    kbold = lin_kernel2(input, Xin, thetad(1:2))' + kernel2(input, Xin, thetad(3:5))';;
    
    A = Xout'*invKd;
    output = A*kbold;
    output = output';
    output_var = zeros(M, 1);
    for i = 1:M
        output_var(i) = lin_kernel2(input(i,:), input(i,:), thetad(1:2)) + ...
            thetad(5) +1/thetad(end) - kbold(:, i)'*invKd*kbold(:, i); 
    end
elseif (modelType(3) == 4)
    [output, output_var] = ...
        lin2ManifoldOutputs(input, Xin, Xout, thetad, invKd);
elseif (modelType(3) == 5)
    N = size(Xout, 1);
    M = size(input, 1);
    
    alpha = zeros(N, q);
    
    kbold = lin_kernel(input, Xin, thetad(1))' + kernel(input, Xin, thetad(2:3))';;
    
    A = Xout'*invKd;
    output = A*kbold;
    output = output';
    output_var = zeros(M, 1);
    for i = 1:M
        output_var(i) = lin_kernel(input(i,:), input(i,:), thetad(1)) + ...
            thetad(3) +1/thetad(end) - kbold(:, i)'*invKd*kbold(:, i); 
    end
end

if (modelType(2) == 0)
    xo = output;
    varxo = output_var;
elseif (modelType(2) == 1)
    xo = output + xi(1:q);
    varxo = output_var;
end
