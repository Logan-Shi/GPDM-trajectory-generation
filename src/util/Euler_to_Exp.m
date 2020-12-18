function v = Euler_to_Exp(x) 

Rx = AxisAngle_2_Rotation([1 0 0 (x(1)/180)*pi]);
Ry = AxisAngle_2_Rotation([0 1 0 (x(2)/180)*pi]);
Rz = AxisAngle_2_Rotation([0 0 1 (x(3)/180)*pi]);

R = Rx*Ry*Rz;
if (x(1) < 1e-6 && x(2) < 1e-6 && x(3) < 1e-6)
    v = zeros(3,1);
elseif (x(1) < 1e-6 && x(2) < 1e-6)
    v = [0; 0; (x(3)/180)*pi];
elseif (x(1) < 1e-6 && x(3) < 1e-6)
    v = [0; (x(2)/180)*pi; 0];
elseif (x(2) < 1e-6 && x(3) < 1e-6)
    v = [(x(1)/180)*pi; 0; 0];
else
    rout = Rotation_2_AxisAngle(R);
    v = rout(1:3);
    v = (v/norm(v))*rout(4);
    v = v';
end