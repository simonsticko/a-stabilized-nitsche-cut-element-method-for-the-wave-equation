function[folderName]=sliverTimeStep()
addpath('assembling')
cclock=clock();
folderName=[num2str(cclock(1)) '_' num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
resultdir=['results/timestep/' folderName '/'];
mkdir(resultdir);
n=flipud(2.^(2:2:6)');
epsilon=2.^(-15:1:-1);
h=zeros(length(n),1);
eigNoStab=zeros(length(n),length(epsilon));
eigOnlyM=zeros(length(n),length(epsilon));
eigOnlyA=zeros(length(n),length(epsilon));
eigStab=zeros(length(n),length(epsilon));
for i=1:length(n)
    for j=1:length(epsilon)
            disp(['i=' num2str(i) ' j=' num2str(j) ', time=']);
            disp(num2str(clock()));
        [h(i),eigNoStab(i,j),eigOnlyM(i,j),eigOnlyA(i,j),eigStab(i,j)]=solveOnce(n(i),epsilon(j));
    end
end
save([resultdir 'savedData']);
end



function[h,eigNoStab,eigOnlyM,eigOnlyA,eigStab]=solveOnce(n,epsilon)
xLim=[-1.5,1.5];
yLim=xLim;
H=diff(xLim)/(n-1);
xPolygon=[xLim(2)-(2-epsilon)*H;2];
yPolygon=[-2;2];
[xB,yB]=meshgrid(xPolygon,yPolygon);
order=convhull(xB,yB);
XB=[xB(:),yB(:)];
XB=XB(order,:);
haveInnerProblem=false;
cutMesh=CutMesh(xLim,yLim,n,XB,haveInnerProblem);
h=cutMesh.h;
gD=@(x,y) zeros(size(x));
f=@(x,y) zeros(size(x));
dirichletInner=true;
dirichletOuter=true;
%Need to choose
[M,A,~,~,~,m,a]=assemble(cutMesh,f,gD,dirichletInner,gD,dirichletOuter);
eigNoStab=eigs(a,m,1);
eigOnlyM=eigs(a,M,1);
eigOnlyA=eigs(A,m,1);
eigStab=eigs(A,M,1);
end
