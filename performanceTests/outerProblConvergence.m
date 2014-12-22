%Makes a convergence test by for the plane wave outer problem by 
%refining the grid, plots relevant quantities and saves them in a folder. 
function[]=outerProblConvergence()
close all;
%make folder to put results in.
cclock=clock();
folderName=[num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
path=['results/outerProbl/' folderName '/'];
mkdir(path);
%n is "fineness" in space want to start with the most computational heavy
%simulation therefore flipud
n=flipud(10*(7:1:12)');
Tend=1;
%Return error, want to save them.
h=zeros(size(n));
bcError=zeros(size(n));
uError=zeros(size(n));
gradError=zeros(size(n));
dudtError=zeros(size(n));
for j=1:length(n)
    disp(['started with j=' num2str(j) ', time=']);
    disp(num2str(clock()));
    [h(j),uError(j),gradError(j),dudtError(j),bcError(j)]=...
        outerWave(n(j),Tend,max(n),false);
end
save([path 'savedData' '.mat']);
%error in u
yLab='$\left\Vert  u_{h}-u \right\Vert $';
plotAndSave(h,uError,yLab,path);
%error in grad
yLab='$\left\Vert \nabla( u_{h}-u )\right\Vert$';
plotAndSave(h,gradError,yLab,path);
%dudt
yLab='$\left\Vert \dot{u}_{h}-\dot{u}\right\Vert$';
plotAndSave(h,dudtError,yLab,path);
%boundary error
yLab='$ \left \Vert \hat{n} \cdot \nabla u_{h} - g_{N} \right \Vert_{\Gamma_{N}}$';
plotAndSave(h,bcError,yLab,path);
end

%Resizes and saves figure.
function[]=plotAndSave(h,error,yLab,path)
fig=figure();
name=inputname(2);
pString='$p$';
color='b';
plotLogarithmic(h,error,'h',yLab,pString,color);
resizeFig();
fixPaperSize();
saveas(figure(fig),[path name '.fig'],'fig');
saveas(figure(fig),[path name '.pdf'],'pdf');
end
