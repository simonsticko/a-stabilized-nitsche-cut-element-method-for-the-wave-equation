function[]=halfperiodic()
addpath analyticFunctions/plane;
addpath assembling/;
n=20;
epsilon=.5;
waveNumber=2*pi/3;
xLim=[-1.5,1.5];
yLim=xLim;
H=diff(xLim)/(n-1);
xPolygon=[-2;2];
yPolygon=[yLim(1)+epsilon*H;yLim(2)-epsilon*H];
[xB,yB]=meshgrid(xPolygon,yPolygon);
order=convhull(xB,yB);
XB=[xB(:),yB(:)];
XB=XB(order,:);
haveInnerProblem=true;
cutMesh=CutMesh(xLim,yLim,n,XB,haveInnerProblem);
h=cutMesh.h;
nodes=create_sets(cutMesh.dt,h);
c=1;
uConst=0;
[M,A,L,u0,du0dt]=createSystem(cutMesh,nodes,c,waveNumber);
cfl=0.1;
domainWidth=diff(xLim);
nMax=n;
periodTime=2*pi/(waveNumber*c);
nPeriods=2;
T=periodTime*nPeriods;
[dt,nSteps]=getnSteps(cfl,c,T,nMax,domainWidth);
[u,dudt,t]=timeStepRK(u0,du0dt,A,M,L,nSteps,dt,c);
u=convertToFull(u,nodes);
dudt=convertToFull(dudt,nodes);
saveName=['periodic_' '.mat'];
save(saveName,'t','u','dudt','cutMesh','xLim','yLim','c',...
    'dt','uConst');
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
[u0,du0dt]=setupInitialConditions(cutMesh,AGlobal,dirichletInner,...
    dirichletOuter,nodes,waveSpeed,waveNumber);
L=@(t) zeros(size(A,1),1);
end

function[u0,du0dt]=setupInitialConditions(cutMesh,A,dirichletInner,dirichletOuter,nodes,c,k)
t=0;
uAnaly=@(x,y) planeWave(x,y,t,k,c);
% graduAnaly=@(x,y) gradPlaneWave(x,y,t,k,c);
dudtAnaly=@(x,y) dplaneWavedt(x,y,t,k,c);
% graddudtAnaly=@(x,y) dgradPlaneWavedt(x,y,t,k,c);
% Au0=getAuRitz(cutMesh,uAnaly,graduAnaly,...
%     dirichletInner,dirichletOuter);
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

function[nodes]=create_sets(triangulation,h)
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