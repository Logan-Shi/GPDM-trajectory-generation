function[Xt] = initMissing(Xt, missing, Xin, Xout, invKd, thetad, modelType)

M = size(Xt,1);
q = size(Xt,2);

for m=1:M
    if (ismember(m, missing))
        if (m == 1)
            Xt(m,:) = zero(1,q);
        end
        if (m == 2)
            if (modelType(1) < 2)
                Xt(m,:) = Xt(m-1,:);
            else
                [Xt(m,:), Xm_var] = ...
                    priorManifoldOutputs([Xt(m-1,:)], ...
                    Xin, Xout, thetad, invKd, modelType);
            end
 
        end
        if (m > 2)
            [Xt(m,:), Xm_var] = ...
                priorManifoldOutputs([Xt(m-1,:) Xt(m-2,:)], ...
                Xin, Xout, thetad, invKd, modelType);
            
        end
    end
end

