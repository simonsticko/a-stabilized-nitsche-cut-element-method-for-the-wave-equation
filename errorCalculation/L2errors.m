%Integrates in polar coordinates in order to obtain the L2-errors of u dudt
%and grad u. theta is integrated from 0 to 2*pi. r is integrated from Rmin
%to Rmax which may be functions of theta. Other inputs are object of class
%uInterpolator and analytic functions.
function[l2Norm,l2NormGrad,l2Normdudt]=...
    L2errors(uInterpol,uAnaly,graduAnaly,dudtAnaly,Rmin,Rmax)
AbsTol=1E-3;
RelTol=1E-2;
%error between analytical and numerical solution squared
f=@(x,y) (uAnaly(x,y)-uInterpol.evaluate(x,y)).^2;
%Transform to polar coordinates and perform the intetral.
%Note that an r comes from transformation to polar coordinates.
integrand=@(theta,r) r.*f(r.*cos(theta),r.*sin(theta));
l2Norm=sqrt(integral2(integrand,0,2*pi,Rmin,Rmax,'AbsTol',AbsTol,'RelTol',RelTol));
disp('Done u error');
%Gradient error
fgrad=@(x,y) gradDiffSq(x,y,graduAnaly,@(xx,yy) uInterpol.evaluategrad(xx,yy));
gradIntegrand=@(theta,r) r.*fgrad(r.*cos(theta),r.*sin(theta));
l2NormGrad=sqrt(integral2(gradIntegrand,0,2*pi,Rmin,Rmax,...
    'AbsTol',10*AbsTol,'RelTol',.1*RelTol));%used 100 and 0.1 previously
disp('Done grad(u) error');
%dudt error:
fdudt=@(x,y) (dudtAnaly(x,y)-uInterpol.evaluatedudt(x,y)).^2;
dudtIntegrand=@(theta,r) r.*fdudt(r.*cos(theta),r.*sin(theta));
l2Normdudt=sqrt(integral2(dudtIntegrand,0,2*pi,Rmin,Rmax,...
    'AbsTol',AbsTol,'RelTol',RelTol));
disp('Done dudt error');
end

function[diffSq]=gradDiffSq(x,y,gradUAnaly,graduNum)
sizeX=size(x);
x=x(:);
y=y(:);
diffGrad=gradUAnaly(x,y)-graduNum(x,y);
%The 2 means take column wise scalar product.
diffSq=dot(diffGrad,diffGrad,2);
diffSq=reshape(diffSq,sizeX);
end