% Engin Tola
% etola@yahoo.com
function R = Quaternion_2_Rotation(q)

x=q(1);
y=q(2);
z=q(3);
w=q(4);


R(1,1) = 1-2*y*y-2*z*z;
R(2,2) = 1-2*x*x-2*z*z;
R(3,3) = 1-2*x*x-2*y*y;

R(2,1) = 2*x*y-2*w*z;
R(3,1) = 2*x*z+2*w*y;

R(1,2) = 2*x*y+2*w*z;
R(1,3) = 2*x*z-2*w*y;

R(3,2) = 2*y*z-2*w*x;
R(2,3) = 2*y*z+2*w*x;