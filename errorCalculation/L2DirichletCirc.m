%Calculates the L2 error between the solution in uInterpol given
%dirichlet boundary condition gD by taking the line integral around the 
%circle with radius R.
function[error]=L2DirichletCirc(uInterpol,R,gD)
% AbsTol=1E-6;
% RelTol=1E-4;
f=@(x,y) (uInterpol.evaluate(x,y)-gD(x,y)).^2;
integrand=@(theta) R*f(R*cos(theta),R*sin(theta));
error=sqrt(integral(integrand,0,2*pi));%,...
%     'AbsTol',AbsTol,'RelTol',RelTol));
end