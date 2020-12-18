function printHPs(model)

load(model);
[K invK] = computeKernel(X, theta);
[Xin Xout] = priorIO(X, segments, modelType);
[Kd invKd] = computePriorKernel(Xin, thetad, modelType(3));
fprintf('log(beta): %4.8f\n', sum(log(theta)));
fprintf('log(alpha): %4.8f\n', sum(log(thetad)));
% fprintf('D/2*logdet(K): %4.8f\n', D/2*logdet(K));
% L = 0;
% for d= 1:D
%     L = L -0.5*w(d)*w(d)*Y(:, d)'*invK*Y(:, d);
% end
% fprintf('-L reconstruction: %4.8f\n', -L - N*sum(log(w)));
% fprintf('q/2*logdet(Kd): %4.8f\n', q/2*logdet(Kd));
% Lp = 0;
% for d= 1:q
%         Lp = Lp - 0.5*Xout(:, d)'*invKd*Xout(:, d);
% end
%     
%     if (modelType(1) == 0 || modelType(1) == 1)
%         XDiff = X(segments+1,:) - X(segments,:);
%         Lp = Lp - 0.5*sum(sum(XDiff.*XDiff));
%     end
%     
%     Lp = Lp - 0.5*sum(sum(X(segments,:).*X(segments,:)));
%     fprintf('-L dynamics: %4.8f\n', -Lp);
% 
% fprintf('-ln(p(Y|X)): %4.8f\n', -(L-D*N/2*log(2*pi)-D/2*logdet(K)+N*sum(log(w))));
% fprintf('-ln(p(X)): %4.8f\n', -(Lp-q*size(Xin,1)/2*log(2*pi)-q/2*logdet(Kd)));
if (modelType(3) == 0)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('Prior Kernel width1: %4.8f\n', 1/thetad(1))
    fprintf('Prior Kernel width2: %4.8f\n', 1/thetad(2))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(3))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(4))
elseif (modelType(3) == 1)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('Prior Kernel width: %4.8f\n', 1/thetad(1))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(2))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(3))
elseif (modelType(3) == 2)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('Prior Process variance: %4.8f\n', thetad(1))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(2))
elseif (modelType(3) == 3)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('Prior Linear Process variance1: %4.8f\n', thetad(1))
    fprintf('Prior Linear Process variance2: %4.8f\n', thetad(2))
    fprintf('Prior Kernel width1: %4.8f\n', 1/thetad(3))
    fprintf('Prior Kernel width2: %4.8f\n', 1/thetad(4))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(5))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(6))
elseif (modelType(3) == 4)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('Prior Process variance1: %4.8f\n', thetad(1))
    fprintf('Prior Process variance2: %4.8f\n', thetad(2))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(3))
elseif (modelType(3) == 5)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('Prior Linear Process variance: %4.8f\n', thetad(1))
    fprintf('Prior Kernel width: %4.8f\n', 1/thetad(2))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(3))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(4))
end
