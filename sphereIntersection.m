function [ i1, i2, err ] = sphereIntersection( p1, p2, p3, r1, r2, r3 )
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/239659
%SPHEREINTERSECTION Summary of this function goes here
i1 = 0;
i2 = 0;
err = 0;
 p21 = p2-p1;
 p31 = p3-p1;
 c = cross(p21,p31);
 c2 = sum(c.^2);
 u1 = cross(((sum(p21.^2)+r1^2-r2^2)*p31 - ...
             (sum(p31.^2)+r1^2-r3^2)*p21)/2,c)/c2;
 v = sqrt(r1^2-sum(u1.^2))*c/sqrt(c2);
 if isreal(v)
     i1 = p1+u1+v;
     i2 = p1+u1-v;
 else
     err = 1;
 end
end

