%Makes a convergence test by for the vibrating membrane problem by 
%refining the grid, plots relevant quantities and saves them in a folder. 
function[]=membraneConvergence()
close all;
%make folder to put results in.
cclock=clock();
folderName=[num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
path=['results/waveDirichlet/' folderName '/'];
mkdir(path);
%n is "fineness" in space want to start with the most computational heavy
%simulation therefore flipud
n=flipud(10*(9:2:17)');
nPeriods=1.125;
h=zeros(size(n));
uError=zeros(size(n));
condM=zeros(size(n));
boundaryError=zeros(size(n));
XMaxError=zeros(length(n),2);
gradError=zeros(size(n));
dudtError=zeros(size(n));
%timestep and number of steps
% nSteps=max(ceil(1+C*nPeriods*n))*ones(size(n));
for j=1:length(n)
    disp(['started with j=' num2str(j) ', time=']);
    disp(num2str(clock()));
    [uError(j),~,h(j),condM(j),boundaryError(j),XMaxError(j,:),...
        gradError(j),dudtError(j)]=fictitiousWave(nPeriods,n(j),max(n),false);
end
pString='$p$';
color='b';
%Plot errror
figure(1);
hold on;
yLab='$\left\Vert  u_{h}-u \right\Vert $';
plotLogarithmic(h,uError,'h',yLab,pString,color);
fixAndSave(path,'L2Error');
%Plot condition number
figure(2);
hold on;
plot(log(h),condM,[color 'o'],'linewidth',2);
xlabel('log(h)');
ylabel('cond(M)');
fixAndSave(path,'condM');
%Plot boundaryError
figure(3);
hold on;
yLab='$\left\Vert u_{h}-g_{D}\right\Vert _{\Gamma_{D}}$';
plotLogarithmic(h,boundaryError,'h',yLab,pString,color);
fixAndSave(path,'boundaryError');
%Plot where maximum error occurs.
figure(4);
hold on;
rMaxError=sqrt(XMaxError(:,1).^2+XMaxError(:,2).^2);
plot(h,rMaxError,[color 'o'],'linewidth',2);
xlabel('h');
ylabel('r of max error');
fixAndSave(path,'maxError');
%Plot boundaryError
figure(5);
hold on;
yLab='$\left\Vert \nabla( u_{h}-u )\right\Vert$';
plotLogarithmic(h,gradError,'h',yLab,pString,color);
fixAndSave(path,'gradError');
%Plot dudtError
figure(6);
hold on;
yLab='$\left\Vert \dot{u}_{h}-\dot{u}\right\Vert$';
plotLogarithmic(h,dudtError,'h',yLab,pString,color);
fixAndSave(path,'dudtError');
save([path 'savedData' '.mat']);
close all;
end

%Resizes and saves figure.
function[]=fixAndSave(path,name)
resizeFig();
fixPaperSize();
saveas(gcf,[path name '.fig'],'fig');
saveas(gcf,[path name '.pdf'],'pdf');
end

