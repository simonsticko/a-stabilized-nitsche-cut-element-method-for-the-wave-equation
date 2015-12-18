function[]=sliverTimeStep()
addpath('assembling')
resultdir='results/timestep/';
n=[20,40];
epsilon=linspace(.5,.7,5);
h=zeros(length(n),1);
stabilized=zeros(length(n),length(epsilon));
unstabilized=zeros(size(stabilized));
for i=1:length(n)
    for j=1:length(epsilon)
        [h(i),stabilized(i,j),unstabilized(i,j)]=solveOnce(n(i),epsilon(j));
    end
end
save([resultdir 'timestep']);
end

function[h,eigStabilized,eigUnstabilized]=solveOnce(n,epsilon)
xLim=[-1.5,1.5];
yLim=xLim;
h=diff(xLim)/(n-1);
xPolygon=[xLim(2)-(2-epsilon)*h;2];
yPolygon=[-2;2];
[xB,yB]=meshgrid(xPolygon,yPolygon);
order=convhull(xB,yB);
XB=[xB(:),yB(:)];
XB=XB(order,:);
haveInnerProblem=false;
cutMesh=CutMesh(xLim,yLim,n,XB,haveInnerProblem);
gD=@(x,y) zeros(size(x));
f=@(x,y) zeros(size(x));
dirichletInner=true;
dirichletOuter=false;
%Need to choose
[M,A,~,~,~,m,a]=assemble(cutMesh,f,gD,dirichletInner,gD,dirichletOuter);
eigStabilized=eigs(A,M,1);
eigUnstabilized=eigs(a,m,1);
end
