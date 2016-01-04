function[]=domainDepedenceGammaM()
close all;
addpath assembling/;
cclock=clock();
folderName=[num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
folder=['results/domainGammaM/' folderName '/'];
mkdir(folder);
n=2^7;
plotMinOnce(folder,n);
findMinOnAxisRatio(n,folder)
end

function[]=plotMinOnce(folderName,n)
axisRatio=1;
[h,m,J]=assembleOnce(n,axisRatio);
condM=@(gammaM) condest(m+gammaM*h^2*J);
fig=figure();
handle=ezplot(condM,[1e-5, 1]);
set(handle,'linewidth',2);
title('');
xlabel('$\gamma_M$','interpreter','latex');
ylabel('cond$\mathcal(M)$','interpreter','latex');
yLim=ylim();
ylim([yLim(1) 300]);
resizeFig(16);
fixPaperSize();
saveas(fig,[folderName 'gammaMminimumCurve.pdf'],'pdf')
saveas(fig,[folderName 'gammaMminimumCurve.fig'],'fig')
end


function[]=findMinOnAxisRatio(n,folderName)
axisRatio=linspace(.1,1,20);
minGammaM=zeros(size(axisRatio));
for i=1:length(axisRatio)
    disp(['started with i=' num2str(i) ', time=']);
    disp(num2str(clock()));
    [h,m,J]=assembleOnce(n,axisRatio(i));
    condM=@(gammaM) condest(m+gammaM*h^2*J);
    minGammaM(i)=fminbnd(condM,1e-5,1);
end
save([folderName 'saveddata.mat']);
%plot minimum.
fig=figure();
plot(axisRatio,minGammaM,'bo')
xlabel('$a/b$','interpreter','latex');
ylabel('$\arg \min_{\gamma_M} cond(m+\gamma_M h^2 J)$','interpreter','latex');
resizeFig(16);
fixPaperSize();
savename=[folderName 'minGammaM'];
saveas(fig,[savename '.pdf'],'pdf');
saveas(fig,[savename '.fig'],'fig');
end

function[h,m,J]=assembleOnce(n,axisRatio)
xLim=[-1.5,1.5];
yLim=xLim;
%size of majoraxis
a=1;
b=axisRatio*a;
theta=linspace(0,2*pi,500)';
XB=[a*cos(theta),b*sin(theta)];
haveInnerProblem=true;
cutMesh=CutMesh(xLim,yLim,n,XB,haveInnerProblem);
gD=@(x,y) zeros(size(x));
f=@(x,y) zeros(size(x));
dirichletInner=true;
%Need to choose
[~,~,~,~,~,m,~,J]=assemble(cutMesh,f,gD,dirichletInner);
h=cutMesh.h;
end