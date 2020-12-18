% Engin Tola
% etola@yahoo.com
function out = Quaternion_2_AxisAngle(AA)

theta = 2*acos(AA(4));

s2 = sin(theta/2);

nx = AA(1) / s2;
ny = AA(2) / s2;
nz = AA(3) / s2;

out = [nx ny nz theta];
