function my_visualize(X, Y, invK, theta, YLbls, meanData, activeSet, ...
			visualiseFunction, visualiseModify, varargin)

% GPLVMVISUALISE Visualise the manifold.

global visualiseInfo

symbol{1} = 'r+';
symbol{2} = 'bo';
symbol{3} = 'mx';
symbol{4} = 'g*';
symbol{5} = 'ys';
symbol{6} = 'rd';
symbol{7} = 'bV';
symbol{8} = 'm^';
symbol{9} = 'g<';
symbol{10} = 'y>';

activeX = X(activeSet, :);
activeY = Y(activeSet, :);

numData = size(Y, 1);
if nargin < 6
  YLbls = ones(numData, 1);
end

x1 = linspace(min(X(:, 1))*1.1, max(X(:, 1))*1.1, 30);
x2 = linspace(min(X(:, 2))*1.1, max(X(:, 2))*1.1, 30);
[X1, X2] = meshgrid(x1, x2);
XTest = [X1(:), X2(:)];
[testY, testYVar] = manifoldOutputs(XTest, activeX, activeY, theta, invK);

testYPlot = testY;
testYPlot(find(testYVar>prctile(testYVar(:, 1), 25))) =NaN;

D = size(Y,2);
A = Y(activeSet, :)'*invK;
loglikelihoods = zeros(size(XTest,1),1);
for i=1:size(XTest,1)
    x_i = XTest(i,:);
    kbold = kernel(x_i, activeX, theta)';
    f = A*kbold;
    sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 

    if sigma2 < 1e-6;
       sigma2 = 1e-6;
    end

    yHat_i = 0; %Y(i, :)' - f;

    L = -D/2*log(sigma2) -D/2*log(2*pi) - yHat_i'*yHat_i/(2*sigma2)-x_i*x_i'/2;
    L = -L;

    loglikelihoods(i) = L;
end



  
figure(1)
clf
% Create the plot for the data
clf
visualiseInfo.plotAxes = axes('position', [0.05 0.05 0.9 0.9]);
hold on
[c, h] = contourf(X1, X2, -reshape(loglikelihoods, size(X1)), 128); 
shading flat
colormap gray;
colorbar
gplvmtwoDPlot(X, YLbls, symbol);


%dataSet(1) = plot(X(:, 1), X(:, 2), 'r.')
%set(visualiseInfo.dataSet(1), 'MarkerSize', 10)
%hold on

% Set up the X limits and Y limits of the main plot
xLim = [min(X(:, 1)) max(X(:, 1))];
xSpan = xLim(2) - xLim(1);
xLim(1) = xLim(1) - 0.05*xSpan;
xLim(2) = xLim(2) + 0.05*xSpan;
xSpan = xLim(2) - xLim(1);

yLim = [min(X(:, 2)) max(X(:, 2))];
ySpan = yLim(2) - yLim(1);
yLim(1) = yLim(1) - 0.05*ySpan;
yLim(2) = yLim(2) + 0.05*ySpan;
ySpan = yLim(2) - yLim(1);

set(visualiseInfo.plotAxes, 'XLim', xLim)
set(visualiseInfo.plotAxes, 'YLim', yLim)
visualiseInfo.digitAxes = [];
visualiseInfo.digitIndex = [];

% Set the callback function
set(gcf, 'WindowButtonMotionFcn', 'classVisualise(''move'')')

figure(2)
clf
if length(varargin)>0
  set(gcf, 'menubar', 'none')
  xPixels = 115;
  yPixels = 115;
  set(gcf, 'position', [232 572 xPixels yPixels/varargin{1}(1)*varargin{1}(2)])
  imageAxesa = subplot(1, 1, 1);
  xWidth = varargin{1}(1)/xPixels;
  yHeight = varargin{1}(2)/yPixels;
  set(imageAxesa, 'position', [0.5-xWidth/2 0.5-yHeight/2 xWidth yHeight])
else
  imageAxesa =subplot(1, 1, 1);
end
visData = zeros(1,size(Y, 2));
visData(1) = min(min(Y));
visData(end) = max(max(Y));
visHandle = feval(visualiseFunction, visData, varargin{:});
set(visHandle, 'erasemode', 'xor')
colorMap gray
set(imageAxesa, 'visible', 'off')

%imageAxesa = axes('position', [0.95 0.05 0.05 0.05]);
%visualiseInfo.visHandle = feval(visualiseFunction, Y(1, :), varargin{:});
%colorMap gray
%set(imageAxesa, 'visible', 'off')

% Pass the data to visualiseInfo
visualiseInfo.X = X(activeSet, :);
visualiseInfo.Y = Y(activeSet, :);
visualiseInfo.theta = theta;
visualiseInfo.invSigma= invK;
visualiseInfo.meanData = meanData;
visualiseInfo.varargin = varargin;
visualiseInfo.visualiseModify = visualiseModify;
visualiseInfo.visHandle = visHandle;
hold off




function returnVal = gplvmtwoDPlot(X, label, symbol)

% GPLVMTWODPLOT Helper function for plotting the labels in 2-D.

returnVal = [];

if ~isempty(label)
  for i = 1:size(X, 1)
    labelNo = find(label(i, :));
    try 
      returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
    catch
      if strcmp(lasterr, 'Index exceeds matrix dimensions.')
	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
      end
    end
  end
else
  returnVal = plot(X(:, 1), X(:, 2), 'rx');
end