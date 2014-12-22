%integrates over a circle with radius R to obtain the L2-errors.
function[l2Norm,boundaryError,l2NormGrad,l2Normdudt]=...
    L2Norm(uInterpol,R,uAnaly,graduAnaly,dudtAnaly)
AbsTol=1E-3;
RelTol=1E-2;
%error between analytical and numerical solution squared
f=@(x,y) (uAnaly(x,y)-uInterpol.evaluate(x,y)).^2;
%Transform to polar coordinates and perform the intetral.
%Note that an r comes from transformation to polar coordinates.
integrand=@(r,theta) r.*f(r.*cos(theta),r.*sin(theta));
l2Norm=sqrt(integral2(integrand,0,R,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol));
%Test dirichlet BC:
integrandBoundary=@(theta) integrand(R,theta);
boundaryError=sqrt(integral(integrandBoundary,0,2*pi,...
    'AbsTol',0.01*AbsTol,'RelTol',RelTol));
%Gradient error
fgrad=@(x,y) gradDiffSq(x,y,graduAnaly,@(xx,yy) uInterpol.evaluategrad(xx,yy));
gradIntegrand=@(r,theta) r.*fgrad(r.*cos(theta),r.*sin(theta));
l2NormGrad=sqrt(integral2(gradIntegrand,0,R,0,2*pi,...
    'AbsTol',100*AbsTol,'RelTol',RelTol));
if(nargin>4)
    %dudt error:
    fdudt=@(x,y) (dudtAnaly(x,y)-uInterpol.evaluatedudt(x,y)).^2;
    dudtIntegrand=@(r,theta) r.*fdudt(r.*cos(theta),r.*sin(theta));
    l2Normdudt=sqrt(integral2(dudtIntegrand,0,R,0,2*pi,...
        'AbsTol',AbsTol,'RelTol',RelTol));
end
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
