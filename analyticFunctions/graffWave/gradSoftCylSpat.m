%Returns spatial part of the analytical solution of the outer problem.
%The full time dependent solution will be.
% u =@(t) real(F*exp(1i*gamma*t));
% dudt =@(t) real(1i*gamma*F*exp(1i*gamma*t));
%This problem is analytically correct up to machine precision
%if gamma = 4*pi if using 145 terms;
function [gradF] = gradSoftCylSpat(X,Y,gamma)
[r,theta,nTerms]=softCommon(X,Y);
R=1;
dFdtheta=0;
%F0 has no dependece on theta.
dFdr=get_dRndr(0,r,gamma,R);
for n = 1:nTerms
    dFndr=get_dRndr(n,r,gamma,R).*getThetan(n,theta);
    dFndtheta=getRn(n,r,gamma,R).*get_dThetandtheta(n,theta);
    dFdtheta = dFdtheta+dFndtheta;
    dFdr = dFdr+dFndr;
end
%see beta page 251.
s=sin(theta);
c=cos(theta);
gradF=[c.*dFdr, s.*dFdr]+[-s.*dFdtheta./r, c.*dFdtheta./r];
end

function[dRndr]=get_dRndr(n,r,gamma,R)
[an]=get_ansoft(n,gamma,R);
dRndr=gamma*(dJdr(n,gamma*r)+an*dHdr(n,gamma*r));
end

function[Thetan]=get_dThetandtheta(n,theta)
Thetan=-n*2*(-1i)^n*sin(n*theta);
end
