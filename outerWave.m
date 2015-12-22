%Solves the problem of a plane wave coming in towards a circle with
%homogenous neumann bc.
%inparamters
%n-mesh parameter
%T end time,
%nMax, maximum n when doing the convergence study.
%ssave boolean defining if the results should be saved in a separate
%mat-files.
%out paramters are
%h mesh size
%L2 errors in u, gradient of u, dudt and 1D-L2-error of the neumann
%boundary condition.
function[h,uError,gradError,dudtError,bcError]=outerWave(n,T,nMax,ssave)
%width of domain.
if(0==nargin)
    n=21;
    T=5;
    nMax=n;
    ssave=true;
end
addpath assembling/;
addpath analyticFunctions/planeWave;
addpath errorCalculation/;
dirichletInner=false;
dirichletOuter=true;
%forcing function for this problem: d2udt2=c^2*(d2udx2+d2udy2)+f(x,y)
f=@(x,y) zeros(size(x));
%Wave speed.
c=1;
%Radius of immersed disc.
R=1;
%Width of domain.
Dwidth=2.11*R;
cfl=0.6;
[dt,nSteps]=getnSteps(cfl,c,T,nMax,2*Dwidth);
%parameter for analytic solution.
gamma=pi;
%Constant shift upwards:
uConst=0;
%Functions for boundary conditions:
gInner=@(x,y) zeros(size(x));
gOutcSq=@(x,y) c^2*softCylSpat(x,y,gamma);
%time dependence
Tdep=@(t) exp(1i*gamma*t);
%Analytic solutions, needed for ritz-projection and for checking the
%l2-errors at the end of the calculation.
uAnaly=@(x,y,t) real(softCylSpat(x,y,gamma)*Tdep(t));
graduAnaly=@(x,y,t) real(gradSoftCylSpat(x,y,gamma)*Tdep(t));
dudtAnaly=@(x,y,t) real(1i*gamma*softCylSpat(x,y,gamma)*Tdep(t));
graddudtAnaly=@(x,y,t) real(1i*gamma*gradSoftCylSpat(x,y,gamma)*Tdep(t));
%Gradient ofanalytic solution of the problem at time equal to zero.
[XB]=getBoundaryPolygon(R);
xLim=[-Dwidth,Dwidth];
yLim=xLim;
cutMesh=CutMesh(xLim,yLim,n,XB,false);
disp('Mesh generated');
% [dt,nSteps]=getnSteps(cfl,cutMesh.h,c,T);
%assemble, LOutX indicates that this is the spatial part of Lout.
[M,A,L,LOutX]=assemble(cutMesh,f,gInner,dirichletInner,gOutcSq,dirichletOuter);
disp('System assembled');
%Complete L
L=@(t) L+real(LOutX*Tdep(t));
%Take Ritz-projection to get initial conditions.
Au0=getAuRitz(cutMesh,@(x,y) uAnaly(x,y,0),@(x,y) graduAnaly(x,y,0),...
    dirichletInner,dirichletOuter);
u0=A\Au0;
disp('u0 Projected');
Adu0dt=getAuRitz(cutMesh,@(x,y) dudtAnaly(x,y,0),@(x,y) graddudtAnaly(x,y,0),...
    dirichletInner,dirichletOuter);
du0dt=A\Adu0dt;
disp('du0dt Projected');
%solve;
[u,dudt,t]=timeStepRK(u0,du0dt,A,M,L,nSteps,dt,c);
disp('Done with timestepping');
%Save solution if required, can use this to play.
if(ssave)
    saveName=['outerWave_' '.mat'];
    save(saveName,'t','u','dudt','cutMesh','xLim','yLim','c',...
        'dt','uConst');
else
    %Maximum integration limit, corresponds to the outer square.
    Rmax=@(theta) Dwidth*min(1./abs(cos(theta)),1./abs(sin(theta)));
    uInterpol=uInterpolator(cutMesh.dt,u(:,end),dudt(:,end),cutMesh.relevant,uConst);
    disp('Created interpolator');
    [bcError]=L2NeumannCirc(uInterpol,R,gInner);
    disp('Calculated Boundary error.');
    [uError,gradError,dudtError]=...
        L2errors(uInterpol,@(x,y) uAnaly(x,y,T),...
    @(x,y) graduAnaly(x,y,T),...
    @(x,y) dudtAnaly(x,y,T),...
    R,Rmax);
end
h=cutMesh.h;
end