%Tests the convergence rate of the ritz-projection.
function[]=RitzConvergence()
addpath assembling/;
addpath analyticFunctions/;
folderName='RitzConv_uConst0';
path=['results/waveDirichlet/' folderName '/'];
mkdir(path);
n=flipud(10*(4:8)');
l2Norm=zeros(size(n));
l2Boundary=zeros(size(n));
l2Grad=zeros(size(n));
h=zeros(size(n));
for j=1:size(n)
    [h(j),l2Norm(j),l2Boundary(j),l2Grad(j)]=solveOnce(n(j));
end
save([path 'savedData']);
savePlot(h,l2Norm,path,'l2error');
savePlot(h,l2Boundary,path,'l2Boundary')
savePlot(h,l2Grad,path,'l2Graderror')
close all;
end

function[]=savePlot(h,y,path,name)
pString='$p$';
color='b';
%Plot errror
fig=figure();
plotLogarithmic(h,y,'h',name,pString,color);
resizeFig;
fixPaperSize;
saveas(fig,[path name '.pdf'],'pdf');
end

function[h,l2Norm,l2Boundary,l2Grad]=solveOnce(n)
%forcing function for this problem: d2udt2=c^2*(d2udx2+d2udy2)+f(x,y)
f=@(x,y) zeros(size(x));
%Radius of domain.
R=1;
%Wave speed.
c=1;
%Which what zero of J0 we want to use.
besselZero=5;
t0=0;
%Constant shift upwards:
uConst=1;
%Functions for boundary conditions:
gD=@(x,y) uConst*ones(size(x));
%Analytic solution of the problem at time equal to zero. 
u0Analy=@(x,y) besselWave(x,y,t0,R,c,besselZero,uConst);
%Gradient ofanalytic solution of the problem at time equal to zero.
gradu0Analy=@(x,y) gradBesselWave(x,y,t0,R,c,besselZero);
%Assemle
[XB]=getBoundaryPolygon(R);
cutMesh=CutMesh(R*[-1.1,1.1],R*[-1.1,1.1],n,n,XB);
[~,A,~]=assemble(cutMesh,f,gD,true);
%Take Ritz-projection to get initial condition.
Au0=getAuRitz(cutMesh,u0Analy,gradu0Analy);
u0=A\Au0;
%Need to now the period of the problem to determine the timestep.
uInterpol=uInterpolator(cutMesh.dt,u0,zeros(size(u0)),cutMesh.relevant,uConst);
[l2Norm,l2Boundary,l2Grad]=...
    L2Norm(uInterpol,R,u0Analy,gradu0Analy);
h=cutMesh.h;
end