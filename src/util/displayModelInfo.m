function displayModelInfo(model)
global USE_GAMMA_PRIOR  % gamma prior for dynamics, only works with RBF kernel
global GAMMA_ALPHA % defines shape of the gamma prior
global USE_LAWRENCE % fix dynamics HPs, as Lawrence suggested (use with thetad = [0.2 0.01 1e6];) 
global FIX_HP % fix all HPs
global MARGINAL_W % marginalize over W while learning X
global MARGINAL_DW % marginalize over scale in dynamics while learning X
global LEARN_SCALE % use different scales for different output dimensions
global REMOVE_REDUNDANT_SCALE % let W absorb the overall scale of reconstruction
global W_VARIANCE % kappa^2 in the paper, not really the variance though
global M_CONST % M value in Jack's master's thesis
global BALANCE % Constant in front of dynamics term, set to D/q for the B-GPDM
global SUBSET_SIZE % Number of data to select for EM, set -1 for all data. 

M_CONST = 1; 
REMOVE_REDUNDANT_SCALE = 1;
LEARN_SCALE = 1; 
MARGINAL_W = 0; 
MARGINAL_DW = 0; 
W_VARIANCE = 1e6; 
FIX_HP = 0; 
USE_GAMMA_PRIOR = 0; 
GAMMA_ALPHA = [5 10 2.5]; 
USE_LAWRENCE = 0;
BALANCE = 1;
SUBSET_SIZE = -1; 

load(model); 
params = [X(:)' log(theta) log(thetad)];
L = -gpdmlikelihood(params, Y, w, segments, modelType); 
fprintf('ln p(Y, X, w, theta, thetad): %4.8f\n', L); 
fprintf('ln p(Y, X | w, theta, thetad): %4.8f\n', L + sum(log(theta)) + sum(log(thetad)) - D*log(2) + D/2*log(2*pi*W_VARIANCE) + 0.5/W_VARIANCE*sum(w.*w)); 
if (modelType(3) == 0)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Kernel width1: %4.8f\n', 1/thetad(1))
    fprintf('Prior Kernel width2: %4.8f\n', 1/thetad(2))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(3))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(4))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt(thetad(3)*thetad(4)));
    fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/(thetad(1) +thetad(2))));
elseif (modelType(3) == 1)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Kernel width: %4.8f\n', 1/thetad(1))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(2))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(3))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt(thetad(2)*thetad(3)));
    fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/thetad(1)));
elseif (modelType(3) == 2)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Process variance: %4.8f\n', thetad(1))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(2))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt(thetad(1)*thetad(2)));
elseif (modelType(3) == 3)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Linear Process variance1: %4.8f\n', thetad(1))
    fprintf('Prior Linear Process variance2: %4.8f\n', thetad(2))
    fprintf('Prior Kernel width1: %4.8f\n', 1/thetad(3))
    fprintf('Prior Kernel width2: %4.8f\n', 1/thetad(4))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(5))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(6))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt((thetad(1) + thetad(2) + thetad(5))*thetad(6)));
    fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/(thetad(3) + theta(4))));
elseif (modelType(3) == 4)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Process variance1: %4.8f\n', thetad(1))
    fprintf('Prior Process variance2: %4.8f\n', thetad(2))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(3))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt((thetad(1)+thetad(2))*thetad(3)));
elseif (modelType(3) == 5)
    fprintf('Kernel width: %4.8f\n', 1/theta(1))
    fprintf('RBF Process variance: %4.8f\n', theta(2))
    fprintf('Noise variance: %4.8f\n', 1/theta(3))
    fprintf('signal to noise ratio: %4.8f\n', sqrt(theta(2)*theta(3)));
    fprintf('characteristic length scale: %4.8f\n', sqrt(1/theta(1)));
    fprintf('Prior Linear Process variance: %4.8f\n', thetad(1))
    fprintf('Prior Kernel width: %4.8f\n', 1/thetad(2))
    fprintf('Prior RBF Process variance: %4.8f\n', thetad(3))
    fprintf('Prior Noise variance: %4.8f\n', 1/thetad(4))
    fprintf('Prior signal to noise ratio: %4.8f\n', sqrt((thetad(1) + thetad(3))*thetad(4)));
    fprintf('Prior characteristic length scale: %4.8f\n', sqrt(1/thetad(2)));
end