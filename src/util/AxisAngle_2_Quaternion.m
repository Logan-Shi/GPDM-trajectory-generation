% Engin Tola
% etola@yahoo.com
function out = AxisAngle_2_Quaternion(inp);

out = [ inp(1)*sin(inp(4)/2) inp(2)*sin(inp(4)/2) inp(3)*sin(inp(4)/2) cos(inp(4)/2)  ];

