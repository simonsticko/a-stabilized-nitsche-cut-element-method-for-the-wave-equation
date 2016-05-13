function[]=convergencePeriodicProblem()
close all;
%make folder to put results in.
cclock=clock();
folderName=[num2str(cclock(1)) '_' num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
path=['results/periodicProblem/convergence' folderName '/'];
mkdir(path);
%n is "fineness" in space want to start with the most computational heavy
%simulation therefore flipud
n=flipud(2.^(3:4)');
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
        calculateErrorsForRefinement(n(j),max(n));
end
save([path 'savedData' '.mat']);
%error in u
yLab='$\left\Vert  u_{h}-u \right\Vert_{\Omega} $';
expectedSlope=2;
plotAndSave(h,uError,yLab,path,expectedSlope);
%error in grad
yLab='$\left\Vert \nabla( u_{h}-u )\right\Vert_{\Omega}$';
expectedSlope=1;
plotAndSave(h,gradError,yLab,path,expectedSlope);
%dudt
yLab='$\left\Vert \dot{u}_{h}-\dot{u}\right\Vert_{\Omega}$';
expectedSlope=2;
plotAndSave(h,dudtError,yLab,path,expectedSlope);
%boundary error
yLab='$ \left \Vert \hat{n} \cdot \nabla u_{h} - g_{N} \right \Vert_{\Gamma_{N}}$';
expectedSlope=1;
plotAndSave(h,bcError,yLab,path,expectedSlope);
end


function[h,uError,gradError,dudtError,boundaryError]=calculateErrorsForRefinement(n,nMax)
saveSolution=false;
waveNumber=2*pi/3;
waveSpeed=1;
nPeriods=1;
[t,u,dudt,cutMesh,xLim,yLim]=halfperiodic(n,nMax,waveSpeed,waveNumber,nPeriods,saveSolution);
endTime=max(t);
uAnaly=@(x,y,t) planeWave(x,y,t,waveNumber,waveSpeed);
dudtAnaly=@(x,y,t) dplaneWavedt(x,y,t,waveNumber,waveSpeed);
uAnalyticEnd=@(x,y) uAnaly(x,y,endTime);
dudtAnalyticEnd=@(x,y) dudtAnaly(x,y,endTime);
graduAnalyEnd=@(x,y) gradPlaneWave(x,y,endTime,waveNumber,waveSpeed);
[uError,gradError,dudtError,boundaryError]=calculateErrors(...
    cutMesh,u(:,end),dudt(:,end),xLim,yLim,...
    uAnalyticEnd,dudtAnalyticEnd,graduAnalyEnd);
h=cutMesh.h;
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
