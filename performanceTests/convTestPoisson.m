%Makes a convergence test by for the vibrating membrane problem by 
%refining the grid, plots relevant quantities and saves them in a folder. 
function[]=convTestPoisson()
close all;
%make folder to put results in.
cclock=clock();
folderName=[num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
path=['../results/poisson/' folderName '/'];
mkdir(path);
%n is "fineness" in space want to start with the most computational heavy
%simulation therefore flipud
n=flipud(10*(4:8)');
%Return errors because I want to save them.
[error,h,boundaryError,XMaxError,gradError]=calculateError(n,'b');
save([path 'savedData' '.mat']);
%Resize and save figures
fixAndSave(path,'L2Error',1);
fixAndSave(path,'boundaryError',3);
fixAndSave(path,'maxError',4);
fixAndSave(path,'gradError',5);
end

%Resizes and saves figure.
function[]=fixAndSave(path,name,fig)
figure(fig);
resizeFig();
fixPaperSize();
saveas(figure(fig),[path name '.fig'],'fig');
saveas(figure(fig),[path name '.pdf'],'pdf');
end

%Loops over the different meshsizes and calculates L2-error on boundary and
%on domain together with conditionnumber and coordinate of where maximum
%error occurs.
function[error,h,boundaryError,XMaxError,gradError]=...
    calculateError(n,color)
h=zeros(size(n));
error=zeros(size(n));
boundaryError=zeros(size(n));
XMaxError=zeros(length(n),2);
gradError=zeros(size(n));
for j=1:length(n)
    disp(['started with j=' num2str(j) ', time=']);
    disp(num2str(clock()));
    [error(j),h(j),boundaryError(j),XMaxError(j,:),...
        gradError(j)]=fictitiousDirichlet(n(j),false);
end
pString='$p$';
%Plot errror
figure(1);
hold on;
yLab='$\left\Vert  u_{h}-u \right\Vert $';
plotLogarithmic(h,error,'h',yLab,pString,color);
%Plot boundaryError
figure(3);
hold on;
plotLogarithmic(h,boundaryError,'h','boundary error',pString,color);
%Plot where maximum error occurs.
figure(4);
hold on;
rMaxError=sqrt(XMaxError(:,1).^2+XMaxError(:,2).^2);
plot(h,rMaxError,[color 'o'],'linewidth',2);
xlabel('h');
ylabel('r of max error');
%Plot boundaryError
figure(5);
hold on;
yLab='$\left\Vert \nabla( u_{h}-u )\right\Vert$';
plotLogarithmic(h,gradError,'h',yLab,pString,color);
end

