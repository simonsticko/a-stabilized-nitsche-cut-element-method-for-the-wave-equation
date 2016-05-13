function[h,uError,gradError,dudtError,boundaryError]=halfperiodic(n,nMax,saveSolution)
if(nargin<1)
    n=22;
    nMax=n;
    saveSolution=true;
end
addpath analyticFunctions/plane;
addpath assembling/;
addpath errorCalculation;
waveNumber=2*pi/3;
xLim=[-1.5,1.5];
[cutMesh,yLim]=createMesh(xLim,n);
nodes=create_sets(cutMesh);
c=1;
uConst=0;
uAnaly=@(x,y,t) planeWave(x,y,t,waveNumber,c);
dudtAnaly=@(x,y,t) dplaneWavedt(x,y,t,waveNumber,c);
[M,A,L,u0,du0dt]=createSystem(cutMesh,nodes,c,waveNumber);
cfl=0.4;
domainWidth=diff(xLim);
periodTime=2*pi/(waveNumber*c);
nPeriods=1;
endTime=periodTime*nPeriods;
[dt,nSteps]=getnSteps(cfl,c,endTime,nMax,domainWidth);
[u,dudt,t]=timeStepRK(u0,du0dt,A,M,L,nSteps,dt,c);
u=convertToFull(u,nodes);
dudt=convertToFull(dudt,nodes);
if(saveSolution)
    saveName=['periodic_' '.mat'];
    save(saveName,'t','u','dudt','cutMesh','xLim','yLim','c',...
        'dt','uConst');
else
    uAnalyticEnd=@(x,y) uAnaly(x,y,endTime);
    dudtAnalyticEnd=@(x,y) dudtAnaly(x,y,endTime);
    graduAnalyEnd=@(x,y) gradPlaneWave(x,y,endTime,waveNumber,c);
    [uError,gradError,dudtError,boundaryError]=calculateErrors(...
        cutMesh,u(:,end),dudt(:,end),xLim,yLim,...
        uAnalyticEnd,dudtAnalyticEnd,graduAnalyEnd);
    h=cutMesh.h;
end
end

function[cutMesh,yDomainLim]=createMesh(xLim,n)
yMeshLim=xLim;
H=diff(xLim)/(n-1);
epsilon=.5;
xPolygon=[-2;2];
yDomainLim=[yMeshLim(1)+epsilon*H;yMeshLim(2)-epsilon*H];
[xB,yB]=meshgrid(xPolygon,yDomainLim);
order=convhull(xB,yB);
XB=[xB(:),yB(:)];
XB=XB(order,:);
haveInnerProblem=true;
cutMesh=CutMesh(xLim,yMeshLim,n,XB,haveInnerProblem);
end

function[M,A,L,u0,du0dt]=createSystem(cutMesh,nodes,waveSpeed,waveNumber)
gD=@(x,y) zeros(size(x));
f=@(x,y) zeros(size(x));
dirichletInner=false;
dirichletOuter=false;%Not used only for an outer problem.
%Need to choose
[MGlobal,AGlobal]=assemble(cutMesh,f,gD,dirichletInner,gD,dirichletOuter);
M=makeMatrixPeriodic(MGlobal,nodes);
A=makeMatrixPeriodic(AGlobal,nodes);
[u0,du0dt]=setupInitialConditions(cutMesh,...
    nodes,waveSpeed,waveNumber);
L=@(t) zeros(size(A,1),1);
end

function[u0,du0dt]=setupInitialConditions(cutMesh,nodes,c,k)
t=0;
uAnaly=@(x,y) planeWave(x,y,t,k,c);
dudtAnaly=@(x,y) dplaneWavedt(x,y,t,k,c);
u0=uAnaly(cutMesh.dt.Points(:,1),cutMesh.dt.Points(:,2));
du0dt=dudtAnaly(cutMesh.dt.Points(:,1),cutMesh.dt.Points(:,2));
u0=removeRightBoundaryDofs(u0,nodes);
du0dt=removeRightBoundaryDofs(du0dt,nodes);
end

function[uFull]=convertToFull(u,nodes)
nInternalNodes=length(nodes.internal);
nBoundaryNodes=length(nodes.leftBoundary);
nNodesInPeriodic=nInternalNodes+nBoundaryNodes;
nNodesInFull=nInternalNodes+2*nBoundaryNodes;
%Make sure u is the right size.
assert(size(u,1)==nInternalNodes+nBoundaryNodes);
internalIndices=1:nInternalNodes;
boundaryIndices=(nInternalNodes+1):nNodesInPeriodic;
uFull=zeros(nNodesInFull,size(u,2));
uFull(nodes.internal,:)=u(internalIndices,:);
uFull(nodes.leftBoundary,:)=u(boundaryIndices,:);
uFull(nodes.rightBoundary,:)=u(boundaryIndices,:);
end

function[u]=removeRightBoundaryDofs(u,nodes)
nInternalNodes=length(nodes.internal);
nBoundaryNodes=length(nodes.leftBoundary);
assert(length(u)==nInternalNodes+2*nBoundaryNodes);
u=[u(nodes.internal);u(nodes.leftBoundary)];
end

function[matrix]=makeMatrixPeriodic(globalMatrix,nodes)
nInternal=length(nodes.internal);
nBoundary=length(nodes.leftBoundary);
nDofs=nInternal+nBoundary;
matrix=sparse(nDofs,nDofs);
upperIndices=1:nInternal;
lowerIndices=(nInternal+1):nDofs;
matrix(upperIndices,upperIndices)=globalMatrix(nodes.internal,nodes.internal);
matrix(upperIndices,lowerIndices)=globalMatrix(nodes.internal,nodes.leftBoundary)+...
    globalMatrix(nodes.internal,nodes.rightBoundary);
matrix(lowerIndices,upperIndices)=globalMatrix(nodes.leftBoundary,nodes.internal)+...
    globalMatrix(nodes.rightBoundary,nodes.internal);
matrix(lowerIndices,lowerIndices)=...
    globalMatrix(nodes.leftBoundary,nodes.leftBoundary)+...
    globalMatrix(nodes.leftBoundary,nodes.rightBoundary)+...
    globalMatrix(nodes.leftBoundary,nodes.rightBoundary)+...
    globalMatrix(nodes.rightBoundary,nodes.rightBoundary);
end

function[nodes]=create_sets(cutMesh)
h=cutMesh.h;
triangulation=cutMesh.dt;
eps=.1*h;
freeBoundaryEdges=triangulation.freeBoundary;
boundaryNodes=unique(freeBoundaryEdges(:));
boundaryPoints=triangulation.Points(boundaryNodes,:);
%Find nodes on the right boundary
isOnRightBoundary=(max(boundaryPoints(:,1))-eps) < boundaryPoints(:,1);
rightBoundary=boundaryNodes(isOnRightBoundary);
%Find nodes on left boundary
isOnLeftBoundary=boundaryPoints(:,1)< (min(boundaryPoints(:,1))+eps);
leftBoundary=boundaryNodes(isOnLeftBoundary);
%Find all nodes in between left and right boundary.
internal=getInternalNodes(triangulation,leftBoundary,rightBoundary);
nodes=struct('leftBoundary',leftBoundary,'rightBoundary',rightBoundary,...
'internal',internal);
end

function[internalNodes]=getInternalNodes(triangulation,leftBoundary,rightBoundary)
allNodes=1:size(triangulation.Points,1);
toKeep=true(size(triangulation.Points,1),1);
toKeep(rightBoundary)=false;
toKeep(leftBoundary)=false;
internalNodes=allNodes(toKeep);
end

function[uError,gradError,dudtError,boundaryError]=calculateErrors(...
    cutMesh,uEnd,dudtEnd,xLim,yLim,uAnaly,dudtAnaly,graduAnaly)
uConst=0;
uInterpol=uInterpolator(cutMesh.dt,uEnd,dudtEnd,cutMesh.relevant,uConst);
AbsTol=1E-3;
RelTol=1E-2;
%Error in u.
integrand=@(x,y) (uAnaly(x,y)-uInterpol.evaluate(x,y)).^2;
uError=sqrt(integral2(integrand,xLim(1),xLim(2),yLim(1),yLim(2),'AbsTol',AbsTol,'RelTol',RelTol));
%Error in dudt
integranddudt=@(x,y) (dudtAnaly(x,y)-uInterpol.evaluatedudt(x,y)).^2;
dudtError=sqrt(integral2(integranddudt,xLim(1),xLim(2),yLim(1),yLim(2),'AbsTol',AbsTol,'RelTol',RelTol));
%Error in grad
gradIntegrand=@(x,y) gradDiffSq(x,y,graduAnaly,@(xx,yy) uInterpol.evaluategrad(xx,yy));
gradError=sqrt(integral2(gradIntegrand,xLim(1),xLim(2),yLim(1),yLim(2),'AbsTol',AbsTol,'RelTol',RelTol));
%Calculate errors on the boundary
normal=[0;1];
upperLine=@(x) (uInterpol.evaluategrad(x,ones(size(x))*yLim(2))*normal ).^2;
lowerLine=@(x) (uInterpol.evaluategrad(x,ones(size(x))*yLim(1))*normal ).^2;
boundaryIntegrand=@(x) upperLine(x)'+lowerLine(x)';
boundaryError=sqrt(integral(boundaryIntegrand,xLim(1),xLim(2),'AbsTol',AbsTol,'RelTol',RelTol));
end