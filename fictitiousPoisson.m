%Solves the Poisson-problem with the Burman-Hansbo method.
function[l2error,h,boundaryError,XMaxError,gradError]=...
    fictitiousPoisson(n,pplot)
addpath assembling/;
if(0==nargin)
    n=30;
    pplot=true;
end
g0=10;
%RHS-forcing function.
f=@(x,y) 1*ones(size(x));
%Boundary condition.
gD=@(x,y) g0*ones(size(x));
%Radius of circle.
R=1;
dirichlet=true;
[XB]=getBoundaryPolygon(R);
cutMesh=CutMesh(R*[-1.5,1.5],R*[-1.5,1.5],n,XB,true);
[~,A,L]=assemble(cutMesh,f,gD,dirichlet);
xi=A\L;
if(pplot)
    surfWrap(cutMesh,xi,g0);
else
    uAnalytic=@(x,y) uAnalytical(x,y,R,g0);
    uGrad=@(x,y) gradu0Analytic(x,y,R);
    uInterpol=uInterpolator(cutMesh.dt,xi,zeros(size(xi)),cutMesh.relevant,g0);
    %Calculate error in L2Norm
    [l2error,boundaryError,gradError]=L2Norm(uInterpol,R,uAnalytic,uGrad);
    [XMaxError]=getLargestErrorCoordinate(xi,cutMesh.dt,uAnalytic,XB,cutMesh.relevant);
    h=cutMesh.h;
end
end

%Analytical solution for the dirichlet equation
% d^2udx^2+d^2udy^2=-1
% u(r=R)=g0
function[uAnaly]=uAnalytical(x,y,R,g0)
r=sqrt(y.^2+x.^2);
uAnaly=(R^2-r.^2)/4.*(r<=R)+g0;
end

%Analytic gradient.
function[uGrad]=gradu0Analytic(x,y,R)
r=sqrt(y.^2+x.^2);
rHat=[x./r y./r];
absGrad=(r<=R).*-r/2;
uGrad=[rHat(:,1).*absGrad, rHat(:,2).*absGrad];
end