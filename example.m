addpath src/gpdm src/gplvm src/netlab src/util

format long

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
global USE_OLD_MISSING_DATA

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


opt = foptions;
opt(1) = 1;
opt(9) = 0;
if MARGINAL_W == 1
    opt(14) = 100; % total number of iterations
    extItr = 1; 
else
    opt(14) = 10; % rescaling every 10 iterations
    extItr = 100; % do extItr*opt(14) iterations in total
end  

% modelType(1) : input of dynamics
%   0 => [x_t, x_{t-1}]
%   1 => [x_t, x_t - x_{t-1}]
%   2 => [x_t]
% modelType(2) : output of dynamics 
%   0 => x_{t+1} 
%   1 => x_{t+1} - x_t
% modelType(3) : kernel type
%   0 => RBF kernel with weighted dimensions, use with input 0 or 1
%   1 => RBF kernel 
%   2 => Linear kernel
%   3 => weighted Linear kernel + RBF kernel with weighted dimensions, use with
%   input 0 or 1
%   4 => weighted linear kernel
%   5 => linear + RBF

% Learn single walker model from lin+rbf kernel.
modelType = [0 1 5]; 
[Y initY varY segments] = loadMocapData({['07_01.amc']}, [1], [2],[260]);
missing = [];
N = size(Y, 1); D = size(Y, 2);
q = 3; % dimensionality of latent space

% PCA
X = zeros(N, q);
refY = Y; meanData = mean(Y);
Y = Y - repmat(meanData, N, 1);
[v, u] = pca(Y);
v(find(v<0))=0;
X = Y*u(:, 1:q)*diag(1./sqrt(v(1:q)));

% initialize hyperparameters
theta = [1 1 exp(1)];
thetad = [0.9 1 0.1 exp(1)];
w = ones(D,1);

[X theta thetad w] = gpdmfitFull(X, Y, w, segments, theta, thetad, opt, ... 
     extItr, modelType, missing);
save example_model X Y w theta thetad modelType N D q meanData  ...
segments initY varY missing refY;
