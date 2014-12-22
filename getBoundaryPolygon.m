%Returns boundary polygon of a circle with Radius r. n is a parameter
%defining how well the circle should be approximated.
function[XB]=getBoundaryPolygon(r,n)
if(nargin<2)
    n=1000;
end
theta=linspace(0,2*pi,n+1)';
xB=r*cos(theta);
yB=r*sin(theta);
XB=[xB yB];
end