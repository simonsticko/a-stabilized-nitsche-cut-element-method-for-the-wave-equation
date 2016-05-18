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
for i=1:length(waveNumbers)
    [speedCenter(i),speedBoundary(i)]=...
        dispersionForWaveNumber(problem,nPointsInEachDirection,waveNumbers(i));
    disp(['i=' num2str(i)])
end
h=problem.cutMesh.gridPointDistance(1);
omegah=problem.waveSpeed*waveNumbers*h;
plotAndSaveWaveErrors(omegah,speedCenter,speedBoundary,folder);
end

function[]=plotAndSaveWaveErrors(omegah,speedCenter,speedBoundary,folder)
fig=figure();
plot(omegah,speedCenter,'bo',omegah,speedBoundary,'rx');
xlim([0 pi/2]);
fontsize=16;
set(gca(),'FontSize',fontsize)
xlabel('$\xi$','interpreter','latex','FontSize',fontsize+4);
ylabel('$c_\xi$','interpreter','latex','FontSize',fontsize+4);
llegend=legend('$\mathcal{C}_1$','$\mathcal{C}_2$');
llegend.set('interpreter','latex','FontSize',fontsize);
saveas(fig,[folder 'waveSpeeds'],'pdf');
saveas(fig,[folder 'waveSpeeds'],'fig');
end

function[waveNumbers]=getWaveNumbersToTest(problem,nPointsInEachDirection)
highestNumberOfWaves=(nPointsInEachDirection-1)/2;
allWaves=1:highestNumberOfWaves;
nWavesToTest=allWaves(1:6);
intervalLength=diff(problem.xLim);
waveNumbers=2*pi/intervalLength*nWavesToTest;
end

function[centerSpeed,boundarySpeed]=dispersionForWaveNumber(problem,n,waveNumber)
nPeriods=1;
saveSolution=false;
yGridSpacing=problem.cutMesh.gridPointDistance(2);
heightAtCenter=mean(problem.yLim);
heightAtBoundary=problem.yLim(1)+.25*yGridSpacing;
xLim=problem.xLim;
waveSpeed=problem.waveSpeed;
endTime=nPeriods*diff(xLim)/waveSpeed;
[t,u,dudt]=problem.solve(waveNumber,endTime,n,saveSolution);
endTime=t(end);
[startInterpolator,endInterpolator]=getInterpolators(u,dudt,problem.cutMesh);
centerSpeed=speedAtHeight(startInterpolator,endInterpolator,xLim,heightAtCenter,endTime,waveSpeed);
boundarySpeed=speedAtHeight(startInterpolator,endInterpolator,xLim,heightAtBoundary,endTime,waveSpeed);
end

function[startInterpolator,endInterpolator]=...
    getInterpolators(u,dudt,cutMesh)
uConst=0;
startInterpolator=uInterpolator(cutMesh.dt,u(:,1),dudt(:,1),cutMesh.relevant,uConst);
endInterpolator=uInterpolator(cutMesh.dt,u(:,end),dudt(:,end),cutMesh.relevant,uConst);
end

function[numericalSpeed]=speedAtHeight(startInterpolator,endInterpolator,xLim,height,...
endTime,exactSpeed)
nPoints=500;
x=linspace(xLim(1),xLim(end),nPoints);
y=ones(size(x))*height;
f1=startInterpolator.evaluate(x,y);
f2=endInterpolator.evaluate(x,y);
plot(x,f1,x,f2);
xPeriod=diff(xLim);
shift=computeSignalShift(f1,f2,xPeriod);
numericalSpeed=exactSpeed+shift/endTime;
end