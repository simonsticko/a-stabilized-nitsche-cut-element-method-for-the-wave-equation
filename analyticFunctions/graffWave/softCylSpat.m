%Returns spatial part of the analytical solution of the outer problem.
%The full time dependent solution will be.
% u =@(t) real(F*exp(1i*gamma*t));
% dudt =@(t) real(1i*gamma*F*exp(1i*gamma*t));
%This problem is analytically correct up to machine precision
%if gamma = 4*pi if using 145 terms;
function [F] = softCylSpat(X,Y,gamma)
[r,theta,nTerms,ssize]=softCommon(X,Y);
R=1;
%F0 has no dependece on theta.
F=getRn(0,r,gamma,R);
for n = 1:nTerms
    Fn=getRn(n,r,gamma,R).*getThetan(n,theta);
    F = F+Fn;
end
F=reshape(F,ssize);
end