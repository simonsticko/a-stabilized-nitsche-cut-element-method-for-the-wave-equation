%Take a range of frequencies that can be resolved 
%on this grid and compare the dispersion errors at one line at the bottom
%(straight through the first bottom row of elements)
%where the stabilization has an effect and with the dispersion error at the
%line straight through the center of the geomtery. Then plot how the
%speeds depend on frequency.
function[]=dispersionError()
close all;
cclock=clock();
folderName=[num2str(cclock(2)) '_' num2str(cclock(3)) '_' ...
    num2str(cclock(4)) '_' num2str(cclock(5))];
folder=['results/dispersionErrors/' folderName '/'];
mkdir(folder);
%nRefinements should be an odd number.
nPointsInEachDirection=31;
problem=halfperiodic(nPointsInEachDirection);
[waveNumbers]=getWaveNumbersToTest(problem,nPointsInEachDirection);
%Test each frequency and measure the wavespeed at center and at the
%boundary.
speedCenter=zeros(length(waveNumbers),1);
speedBoundary=zeros(size(speedCenter));
amplitudeCenter=zeros(length(waveNumbers),1);
amplitudeBoundary=zeros(size(speedCenter));
%Versions without stabilization
speedCenterNoStab=zeros(length(waveNumbers),1);
speedBoundaryNoStab=zeros(size(speedCenter));
amplitudeCenterNoStab=zeros(length(waveNumbers),1);
amplitudeBoundaryNoStab=zeros(size(speedCenter));
%
solveStabilized=@(waveNumber,endTime,n) problem.solve(waveNumber,endTime,n,false);
solveUnStabilized=@(waveNumber,endTime,n) problem.solveWithoutStabilization(waveNumber,endTime,n,false);
%Take each frequency independently.
for i=1:length(waveNumbers)
    disp(['i=' num2str(i)])
    %Solve with stabilization.
    [speedCenter(i),speedBoundary(i),amplitudeCenter(i),amplitudeBoundary(i)]=...
        dispersionForWaveNumber(problem,nPointsInEachDirection,waveNumbers(i),solveStabilized);
    %Solve without stabilization.
    [speedCenterNoStab(i),speedBoundaryNoStab(i),amplitudeCenterNoStab(i),amplitudeBoundaryNoStab(i)]=...
        dispersionForWaveNumber(problem,nPointsInEachDirection,waveNumbers(i),solveUnStabilized);
end
h=problem.cutMesh.gridPointDistance(1);
omegah=problem.waveSpeed*waveNumbers*h;
plotAndSave(omegah,speedCenter,speedBoundary,speedCenterNoStab,speedBoundaryNoStab,...
'$c_\xi$',folder,'waveSpeed');
plotAndSave(omegah,amplitudeCenter,amplitudeBoundary,amplitudeCenterNoStab,amplitudeBoundaryNoStab,...
'$\frac{A_\xi}{A}$',folder,'amplitudes');
end

function[]=plotAndSave(omegah,speedCenter,speedBoundary,speedCenterNoStab,speedBoundaryNoStab,yLabel,folder,name)
fig=figure();
plot(omegah,speedCenter,'bo',omegah,speedBoundary,'rx',...
    omegah,speedCenterNoStab,'g+',omegah,speedBoundaryNoStab,'k^','linewidth',2);
xlim([0 pi/2]);
fontsize=16;
set(gca(),'FontSize',fontsize)
xlabel('$\xi$','interpreter','latex','FontSize',fontsize+8);
ylabel(yLabel,'interpreter','latex','FontSize',fontsize+8);
llegend=legend('$\mathcal{C}_1$ stabilized','$\mathcal{C}_2$ stabilized',...
    '$\mathcal{C}_1$ unstabilized','$\mathcal{C}_2$ unstabilized');
llegend.set('interpreter','latex','FontSize',fontsize,'location','southWest');
fixPaperSize();
saveas(fig,[folder name],'pdf');
saveas(fig,[folder name],'fig');
end

function[waveNumbers]=getWaveNumbersToTest(problem,nPointsInEachDirection)
highestNumberOfWaves=(nPointsInEachDirection-1)/2;
allWaves=1:highestNumberOfWaves;
nWavesToTest=allWaves(1:6);
intervalLength=diff(problem.xLim);
waveNumbers=2*pi/intervalLength*nWavesToTest;
end

function[centerSpeed,boundarySpeed,centerAmplitude,boundaryAmplitude]=...
    dispersionForWaveNumber(problem,n,waveNumber,solveFunction)
nPeriods=1;
yGridSpacing=problem.cutMesh.gridPointDistance(2);
heightAtCenter=mean(problem.yLim);
heightAtBoundary=problem.yLim(1)+.25*yGridSpacing;
xLim=problem.xLim;
waveSpeed=problem.waveSpeed;
endTime=nPeriods*diff(xLim)/waveSpeed;
[t,u,dudt]=solveFunction(waveNumber,endTime,n);
endTime=t(end);
[startInterpolator,endInterpolator]=getInterpolators(u,dudt,problem.cutMesh);
[centerSpeed,centerAmplitude]=speedAtHeight(startInterpolator,endInterpolator,xLim,heightAtCenter,endTime,waveSpeed);
[boundarySpeed,boundaryAmplitude]=speedAtHeight(startInterpolator,endInterpolator,xLim,heightAtBoundary,endTime,waveSpeed);
end

function[startInterpolator,endInterpolator]=...
    getInterpolators(u,dudt,cutMesh)
uConst=0;
startInterpolator=uInterpolator(cutMesh.dt,u(:,1),dudt(:,1),cutMesh.relevant,uConst);
endInterpolator=uInterpolator(cutMesh.dt,u(:,end),dudt(:,end),cutMesh.relevant,uConst);
end

function[numericalSpeed,amplitudeError]=speedAtHeight(startInterpolator,endInterpolator,xLim,height,...
endTime,exactSpeed)
nPoints=500;
x=linspace(xLim(1),xLim(end),nPoints);
y=ones(size(x))*height;
f1=startInterpolator.evaluate(x,y);
f2=endInterpolator.evaluate(x,y);
% plot(x,f1,x,f2);
xPeriod=diff(xLim);
shift=computeSignalShift(f1,f2,xPeriod);
numericalSpeed=exactSpeed+shift/endTime;
%Since cos(x)^2+sin(x)^2=1, the integral of the f2^2 should be equal to
%domainSize/2. So the root mean square should be equal 1/sqrt(2).
domainSize=diff(xLim);
analyticRMS=1/sqrt(2);
amplitudeError=sqrt(trapz(x,f2.^2)/domainSize)/analyticRMS;
end