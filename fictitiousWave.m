%Solves the wave equation with constant dirichlet bc  on a circular domain 
%with radius R. Starting from a circular symmetric bessel function. 
%That is, solves a vibrating membrane problem.
%Inparameters:
%nPeriods- numer of periods we should calculate.
%n -mesh parameter
%nMax maximum n used in the convergence test.
%ssave a boolean defineing wether or not we should save the results in a
%seperate mat-file.
%Returns L2Error of u timestep, meshsize, mass matrix condition number
%1D-l2error along the boundary, coordinte of maximum error,
%L2-error in the gradient and L2-error of the timederivative
function[uError,dt,h,condM,boundaryError,XMaxError,gradError,dudtError]...
    =fictitiousWave(nPeriods,n,nMax,ssave)
if(0==nargin)
    nPeriods=4.125;
    n=60;
    ssave=true;
    nMax=n;
end
addpath assembling/;
addpath analyticFunctions/besselWave/;
addpath errorCalculation/;
%forcing function for this problem: d2udt2=c^2*(d2udx2+d2udy2)+f(x,y)
f=@(x,y) zeros(size(x));
%Radius of domain.
R=1;
%Wave speed.
c=1;
%Which zero of J0 we want to use.
besselZero=5;
t0=0;
%Constant shift upwards:
uConst=2;
%Functions for boundary conditions:, scale with c^2 
gD=@(x,y) uConst*ones(size(x));
gDcSq=@(x,y) c^2*gD(x,y);
%Analytic solution of the problem at time equal to zero.
u0Analy=@(x,y) besselWave(x,y,t0,R,c,besselZero,uConst);
%Gradient ofanalytic solution of the problem at time equal to zero.
gradu0Analy=@(x,y) gradBesselWave(x,y,t0,R,c,besselZero);
%Assemle
[XB]=getBoundaryPolygon(R);
xLim=R*[-1.51,1.51];
yLim=xLim;
cutMesh=CutMesh(xLim,yLim,n,XB,true);
%use dirichlet bc
dirichlet=true;
[M,A,L]=assemble(cutMesh,f,gDcSq,dirichlet);
%Take Ritz-projection to get initial condition.
Au0=getAuRitz(cutMesh,u0Analy,gradu0Analy,dirichlet);
u0=A\Au0;
dudt0=zeros(size(u0));
%L needs to be a function of time.
L=@(t) L;
%Need to now the period of the problem to determine the timestep.
[~,T]=u0Analy(0,0);
cfl=0.05;
[dt,nSteps]=getnSteps(cfl,c,T*nPeriods,nMax,diff(xLim));
%solve;
[u,dudt,t]=timeStepRK(u0,dudt0,A,M,L,nSteps,dt,c);
%Save solution if required, can use this to play.
if(ssave)
    saveName=['uWaveSolution_' num2str(n) '.mat'];
    save(saveName,'t','u','dudt','cutMesh','xLim','yLim','c',...
        'dt','nPeriods','besselZero','uConst');
else
    %calculate the error at the end of the calculation for all nodes inside the
    %domain.
    %Error in L2Norm at end of calculation.
    uAnalyEnd=@(x,y) besselWave(x,y,t(end),R,c,besselZero,uConst);
    dudtAnalyEnd=@(x,y) ddtbesselWave(x,y,t(end),R,c,besselZero);
    graduAnalyEnd=@(x,y) gradBesselWave(x,y,t(end),R,c,besselZero);
    %Create an interpolation object that is defined for all x and y.
    uInterpol=uInterpolator(cutMesh.dt,u(:,end),dudt(:,end),cutMesh.relevant,uConst);
    %     [L2Error,boundaryError,gradError,dudtError]=...
    %         L2Norm(uInterpol,R,uAnalyEnd,graduAnalyEnd,dudtAnalyEnd);
    [uError,gradError,dudtError]=...
        L2errors(uInterpol,uAnalyEnd,graduAnalyEnd,dudtAnalyEnd,0,R);
    [boundaryError]=L2DirichletCirc(uInterpol,R,gD);
    %Want to know where the largest error occurs get this coordinate.
    [XMaxError]=getLargestErrorCoordinate(u(:,end),cutMesh.dt,uAnalyEnd,XB,cutMesh.relevant);
    condM=condest(M);
    h=cutMesh.h;
end
end
