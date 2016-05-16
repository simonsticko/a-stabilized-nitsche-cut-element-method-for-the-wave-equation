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
saveSolution=false;
waveNumber=2*pi/3;
waveSpeed=1;
nPeriods=1;

[t,u,dudt,cutMesh,xLim,yLim]=halfperiodic(n,nMax,waveSpeed,waveNumber,nPeriods,saveSolution);
uInterpol=
end