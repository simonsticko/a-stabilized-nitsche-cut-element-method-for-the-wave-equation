%The dispersion error is typically not dependent on the grid size (the only
%difference is that one can resolve higher frequencies with a finer grid),
%thus we fix a grid. Then take a range of frequencies that can be resolved 
%on this grid and compare the dispersion errors at one line at the bottom
%(straight through the first bottom row of elements)
%where the stabilization has an effect and with the dispersion error at the
%line straight through the center of the geomtery.
%This would give a figure similar to Fig 1.1 in Gustafsson 
%(doi:10.1007/978-3-540-74993-6), with two lines.
function[]=dispersionError()
close all;


nFrequencies=50;
omegah=linspace(0,pi,nFrequencies+2);
omegah=omegah(2:end-1);
n=50;
problem=halfperiodic(n);
h=diff(problem.xLim)/(n-1);
waveNumbersToTest=omegah/h;

speed=zeros(length(waveNumbersToTest),1);
for i=1:length(waveNumbersToTest)
    speed(i)=dispersionForWaveNumber(problem,n,waveNumbersToTest(i));
    disp(['i=' num2str(i)])
end
plot(omegah,speed,'o');
xlim([0 pi]);
ylim([0 problem.waveSpeed]);
end

function[numericalSpeed]=dispersionForWaveNumber(problem,n,waveNumber)
nPeriods=2;
saveSolution=false;
height=mean(problem.yLim);
xLim=problem.xLim;
waveSpeed=problem.waveSpeed;
[t,u,dudt]=problem.solve(waveNumber,nPeriods,n,saveSolution);
endTime=t(end);
[startInterpolator,endInterpolator]=getInterpolators(u,dudt,problem.cutMesh);
numericalSpeed=speedAtHeight(startInterpolator,endInterpolator,xLim,height,endTime,waveSpeed);
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