%Returns spatial part of the analytical solution of the outer problem.
%The full time dependent solution will be.
% u =@(t) real(F*exp(1i*gamma*t));
% dudt =@(t) real(1i*gamma*F*exp(1i*gamma*t));
%This problem is analytically correct up to machine precision
%if gamma = 4*pi if using 145 terms;
function [Fall] = hardCylSpat(X,Y,gamma)
[r,theta,outside,nTerms,ssize]=hardCommon(X,Y,gamma);
%Generate expansion for all coordinates that are outside.
alpha_0=besselj(1,gamma*1)/besselh(1,2,gamma*1);
F = besselj(0,gamma*r)-alpha_0*besselh(0,2,gamma*r);
for n = 1:nTerms
    alpha_n=(besselj(n+1,gamma*1)-besselj(n-1,gamma*1))/...
        (besselh(n+1,2,gamma*1)-besselh(n-1,2,gamma*1));
    F = F+2*(-1i)^n*cos(n*theta).*(besselj(n,gamma*r)-alpha_n*besselh(n,2,gamma*r));
end
%shape back to orignial size:
Fall=zeros(size(outside));
Fall(outside)=F;
Fall=reshape(Fall,ssize);
end
