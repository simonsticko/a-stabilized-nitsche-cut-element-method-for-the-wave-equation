%Returns the gradient of the spatial part of the analytical solution of the
%outer problem. The full time dependent solution will be.
% u =@(t) real(gradF*exp(1i*gamma*t));
%This problem is analytically correct up to machine precision
%if gamma = 4*pi;
function [gradFfull] = gradHardCylSpat(X,Y,gamma)
[r,theta,outside,nTerms]=hardCommon(X,Y,gamma);
%Generate expansion for all coordinates that are outside.
alpha_0=besselj(1,gamma*1)/besselh(1,2,gamma*1);
dF0dr=-2*gamma*(besselj(1,gamma*r)-alpha_0*besselh(1,2,gamma*r));
dFdr=.5*dF0dr;
%F0 has no theta dependence.
dFdtheta=0;
for n = 1:nTerms
    alpha_n=(besselj(n+1,gamma*1)-besselj(n-1,gamma*1))/...
        (besselh(n+1,2,gamma*1)-besselh(n-1,2,gamma*1));
    %term derivative w.r.t theta
    dFndtheta=-2*n*(-1i)^n*sin(n*theta).*(besselj(n,gamma*r)-alpha_n*besselh(n,2,gamma*r));
    %term derivative w.r.t r
    dFndr=gamma*(-1i)^n*cos(n*theta).*(...
        (besselj(n-1,gamma*r)-besselj(n+1,gamma*r))-...
        alpha_n*(besselh(n-1,2,gamma*r)-besselh(n+1,2,gamma*r)));
    
    dFdr = dFdr+dFndr;
    dFdtheta=dFdtheta+dFndtheta;
end
%see beta page 251.
s=sin(theta);
c=cos(theta);
gradF=[c.*dFdr, s.*dFdr]+[-s.*dFdtheta./r, c.*dFdtheta./r];
gradFfull=zeros(size(outside,1),2);
if(~isempty(gradF))
    gradFfull(outside,:)=gradF;
end
end