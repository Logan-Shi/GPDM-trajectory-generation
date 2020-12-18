% Engin Tola
% etola@yahoo.com
function out = Rotation_2_Quaternion(R)

aa = Rotation_2_AxisAngle(R);
out = AxisAngle_2_Quaternion(aa);