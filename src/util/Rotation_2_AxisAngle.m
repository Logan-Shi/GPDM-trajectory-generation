% Engin Tola
% etola@yahoo.com
function rout = Rotation_2_AxisAngle(R)

txy = ( R(2,1) +  R(1,2) ) /2;
txz = ( R(3,1) +  R(1,3) ) /2;
tyz = ( R(3,2) +  R(2,3) ) /2;

sx = R(2,3) - tyz;
sy = R(3,1) - txz;
sz = R(1,2) - txy;

k = txy / ( sx*sy ); 

s = sqrt( (2*k-1)/(k*k) );

theta = asin(s);


minErr = 1e20;
minInd = 1;

th = [ theta -theta ];

for i = 1:2

    x = sx / sin(th(i));
    y = sy / sin(th(i));
    z = sz / sin(th(i));
    
    rr = AxisAngle_2_Rotation([ x,y,z,th(i) ]);
    
    err = sum(sum(abs(rr-R)));

    if err < minErr
        minErr = err;
        minInd = i;
    end
end

x = sx / sin(th(minInd));
y = sy / sin(th(minInd));
z = sz / sin(th(minInd));

rout = [ x y z th(minInd)];