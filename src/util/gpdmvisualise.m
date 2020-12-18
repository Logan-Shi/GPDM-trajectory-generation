function gpdmvisualise(X, Y, segments, invK, theta, invKd, thetad,  modelType, plotType)

% GPLVMVISUALISE Visualise the manifold.
side_max = max(max(X));
side_min = min(min(X));
oldaxis = [side_min side_max side_min side_max side_min side_max]; 

side_min = min(oldaxis); 
side_max = max(oldaxis);


D = size(Y,2);
q = size(X,2);
x = [];
m = 50;
bigm = 50;


% newaxis = oldaxis; 
% for n = 1:3
%     diff = (side_max - side_min);
%     curdim = newaxis([2*n-1 2*n]); 
%     newaxis([2*n-1 2*n]) = ((bigm-1)/(diff))*curdim + (0.5 - (bigm-1)*side_min/diff);
% end

%halfl = 0.5*(side_max - side_min)/(m-1); 
%slope = (m - 0.5 - 0.5)/(side_max - side_min); 
slope = (bigm - 0.5*bigm/m - 0.5*bigm/m)/(side_max - side_min); 
%con = 0.5 - slope*side_min;  
con = 0.5*bigm/m - slope*side_min; 

for n = 1:3
    x{n} = linspace(side_min, side_max, m);
    newaxis([2*n-1 2*n]) = slope*oldaxis([2*n-1 2*n]) + con; 
     x{n}(x{n} < (oldaxis(2*n-1))) = []; 
     x{n}(x{n} > (oldaxis(2*n))) = []; 
end


[X1, X2, X3] = meshgrid(x{1}, x{2}, x{3});
XTest = [X1(:), X2(:), X3(:)];
Z = zeros(size(X1));

if (plotType == 0)
    [testY, testYVar] = manifoldOutputs(XTest, X, Y, theta, invK);
 
    Z = -D/2*log(reshape(testYVar, size(Z,1), size(Z,2), size(Z,3))) - ... 
        reshape(0.5*sum(XTest.*XTest, 2), size(Z,1), size(Z,2), size(Z,3)); 
    Z = Z - D/2*(2*pi) - 1/2*(2*pi); 

    if exist('W', 'var')
        Z = Z + log(det(W));
    end
elseif (plotType == 1)
    [Xin Xout] = priorIO(X, segments, modelType);
    [testY, testYVar] = priorManifoldOutputs([XTest zeros(size(XTest,1),q)] , Xin, Xout, thetad, invKd, modelType);
    Z = -q/2*log(reshape(testYVar,size(Z,1), size(Z,2), size(Z,3))) - ... 
        reshape(0.5*sum(XTest.*XTest, 2), size(Z,1), size(Z,2), size(Z,3));
    Z = Z - q/2*(2*pi) - 1/2*(2*pi); 
elseif (plotType == 2)
    [Xin Xout] = priorIO(X, segments, modelType);
    %[testY1, testYVar1] = manifoldOutputs(XTest, X, Y, theta, invK);
    [testY2 ,testYVar2] = ...
        priorManifoldOutputs([XTest zeros(size(XTest,1),q)] , Xin, Xout, thetad, invKd, modelType);
    [testY3, testYVar1] = manifoldOutputs(testY2, X, Y, theta, invK);
    Z = -D/2*log(reshape(testYVar1, size(Z,1), size(Z,2), size(Z,3))) - ... 
        q/2*log(reshape(testYVar2, size(Z,1), size(Z,2), size(Z,3))) ... 
        - reshape(0.5*sum(XTest.*XTest, 2), size(Z,1), size(Z,2), size(Z,3));
    Z = Z - D/2*(2*pi) - 1/2*(2*pi);
    if exist('W', 'var')
        Z = Z + log(det(W));
    end
end


cut = []; 
for n = 1:3
    %halfl = 0.5*(side_max - side_min)/(bigm-1);
    x{n} = linspace(side_min, side_max, bigm);
    cut{n} = sum(x{n} < (oldaxis(2*n-1))); 
      x{n}(x{n} < (oldaxis(2*n-1))) = []; 
      x{n}(x{n} > (oldaxis(2*n))) = []; 
end


[X1i, X2i, X3i] = meshgrid(x{1}, x{2}, x{3});

Zi = interp3(X1, X2, X3, Z, X1i, X2i, X3i);

%Zi = Z;
%size(Zi)
%  h = vol3d('cdata',Z,'texture','3D', 'Parent', [min(X(:, 1))*1.1 max(X(:, 1))*1.1 ...
%  min(X(:, 2))*1.1 max(X(:, 2))*1.1 ...
%    min(X(:, 3))*1.1 max(X(:, 3))*1.1]);
h = vol3d('cdata',Zi,'texture','3D');
mask = segments-1;
mask(find(mask == 0)) = [];
%SX = X;
% for n = 1:q
%     diff = (side_max - side_min);
%     SX(:,n) = ((bigm-1)/(diff))*SX(:,n) + (0.5 - (bigm-1)*side_min/diff);
% end
% newaxis
%S = SX;
% 
% if size(S,1) > 1
%     dS = S(2:end,:) - S(1:end-1,:);
% 
%     dS(mask,:) = zeros(length(mask),q);
% 
%     dS = [dS; zeros(1,q)];
% 
%     quiver3(S(:,1), S(:,2), S(:,3), dS(:,1), dS(:,2), dS(:,3), 0, 'b');
% end
% hold on;
%plot3(S(:,1), S(:,2), S(:,3), 'bo');

shading flat
%colormap gray;
% axis tight;  daspect([1 1 .4])
alphamap('rampup');
alphamap(.06 .* alphamap);
%view(3);
% Update view since 'texture' = '2D'
hold on;
vol3d(h);

%view(10,10); 
%saxis = axis
%return;
%return;

SX = X;
for n = 1:3
    %slope = (bigm - 0.5*bigm/m - 0.5*bigm/m)/(side_max - side_min); 
    %con = 0.5*bigm/m - slope*side_min; 
    SX(:,n) = slope*SX(:,n) + con - cut{n}; 
end
S = SX;

if size(S,1) > 1
    dS = S(2:end,:) - S(1:end-1,:);
    %size(zeros(length(mask),q))
    dS(mask,:) = zeros(length(mask),q);

    dS = [dS; zeros(1,q)];

    quiver3(S(:,1), S(:,2), S(:,3), dS(:,1), dS(:,2), dS(:,3), 0, 'b');
end

plot3(S(:,1), S(:,2), S(:,3), 'bo');
newaxis(1:2) = newaxis(1:2) - cut{1}; 
newaxis(3:4) = newaxis(3:4) - cut{2}; 
newaxis(5:6) = newaxis(5:6) - cut{3};
axis(newaxis); 
%axis(newaxis); 
%alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')

%axis([side_min side_max side_min side_max side_min side_max side_min side_max]); 
% grid on;
% axis equal;


