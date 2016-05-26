clear;
nPeriods=10;
n=31;
problem=halfperiodic(n);
intervalLength=diff(problem.xLim);
h=intervalLength/(n-1);
highestNumberOfWaves=(n-1)/2;
nWavesInInterval=1:floor(highestNumberOfWaves/4);
nWaves=nWavesInInterval(end);
waveNumber=2*pi/intervalLength*nWaves;
saveSolution=true;
endTime=nPeriods*diff(problem.xLim)/problem.waveSpeed;
[t,u,dudt]=problem.solve(waveNumber,endTime,n,saveSolution);
playWave('periodic_',.5*[1 1 1],inf,true);
% trisurf(problem.cutMesh.dt.ConnectivityList,problem.cutMesh.dt.Points(:,1),problem.cutMesh.dt.Points(:,2),u(:,end))