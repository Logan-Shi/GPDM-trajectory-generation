% Engin Tola
% etola@yahoo.com
function r = AxisAngle_2_Rotation(axisAngle)

x = axisAngle(1);
y = axisAngle(2);
z = axisAngle(3);
theta = axisAngle(4);

s = sin(theta);
c = cos(theta);
t = 1-c;

r(1,1) = t*x*x+c; 
r(2,2) = t*y*y+c; 
r(3,3) = t*z*z+c; 

r(1,2) = t*x*y+s*z;
r(2,1) = t*x*y-s*z;

r(1,3) = t*x*z-s*y;
r(3,1) = t*x*z+s*y;

r(2,3) = t*y*z+s*x;
r(3,2) = t*y*z-s*x;

return;

