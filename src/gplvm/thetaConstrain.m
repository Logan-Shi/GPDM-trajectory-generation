function [theta, gtheta] = thetaConstrain(theta, gtheta)

% THETACONSTRAIN Prevent kernel parameters from getting too big or small.

minTheta = 1e-6;
maxTheta = 1/minTheta;

if any(theta <= minTheta)
  theta(find(theta<=minTheta)) = minTheta;
  if nargin > 1
      %gtheta(find(theta<=minTheta)) = 0;
  end
end
if any(theta>=maxTheta)
  theta(find(theta>=maxTheta)) = maxTheta;
  if nargin > 1
      %gtheta(find(theta>=maxTheta)) = 0; 
  end
end
