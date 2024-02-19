function [tangent,normal]=rotate2princax(a,b,theta);

% theta = angle of rotation, math notation (east == 0, north=90),degrees
% a = eastward component (along x-axis)
% b = northward component (along y-axis)
% tangent = component along major axis of principal ellipse
% normal = component along minor axis of principal ellipse 

tangent=real((a+sqrt(-1)*b).*exp(-(sqrt(-1))*theta*pi/180));

normal=imag((a+sqrt(-1)*b).*exp(-(sqrt(-1))*theta*pi/180));
