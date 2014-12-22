%Assembles the inner problem for different mesh sizes and plots the condition 
%number of the matrix M agains the mesh size h.
function[]=checkCondM()
close all;
addpath assembling/;
cclock=clock();
folderName=[num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
folder=['results/condM/' folderName '/'];
mkdir(folder);
addpath '..'
saveName=[folder 'condM'];
mkdir(folder);
n=flipud(10*(6:2:30)');
condM=zeros(size(n));
h=zeros(size(n));
for j=1:length(n)
    disp(num2str(clock()));
    disp(['n=' num2str(n(j))]);
    [condM(j),h(j)]=getcondM(n(j));
end
save(saveName);
plot(log(h),condM,'bo');
xlabel('log(h)');
ylabel('cond(M)');
resizeFig();
fixPaperSize();
saveas(gcf,saveName,'pdf')
saveas(gcf,saveName,'fig')
end


function[condM,h]=getcondM(n)
f=@(x,y) zeros(size(x));
%Radius of domain.
R=1;
haveInnerProblem=true;
%Constant shift upwards:
uConst=1;
%Functions for boundary conditions:
gD=@(x,y) uConst*ones(size(x));
dirichlet=true;
[XB]=getBoundaryPolygon(R);
cutMesh=CutMesh(R*[-1.1,1.1],R*[-1.1,1.1],n,XB,haveInnerProblem);
[M]=assemble(cutMesh,f,gD,dirichlet);
condM=condest(M);
h=cutMesh.h;
end
