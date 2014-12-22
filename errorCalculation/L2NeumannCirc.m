%Calculates the L2 error betweeh the solution in uInterpol given
%neumannn boundary condition gN by taking the line integral around the 
%circle with radius R.
function[error]=L2NeumannCirc(uInterpol,R,gN)
% AbsTol=1E-6;
% RelTol=1E-4;
f=@(x,y) integrandxy(x,y,uInterpol,gN);
integrand=@(theta) R*f(R*cos(theta),R*sin(theta));
error=sqrt(integral(integrand,0,2*pi));%,...
%     'AbsTol',AbsTol,'RelTol',RelTol));
end

%integral might call the integrand with matrices while the function
%uInterpol.evaluategrad requires vectors. This function is used to
%circument this, by reshaping and shaping back.
function[integrand]=integrandxy(x,y,uInterpol,gN)
sizeX=size(x);
x=x(:);
y=y(:);
grad=uInterpol.evaluategrad(x,y);
r=sqrt(x.^2+y.^2);
normal=[x./r,y./r];
nGrad=sum(grad.*normal,2);
integrand=(nGrad-gN(x,y)).^2;
integrand=reshape(integrand,sizeX);
end