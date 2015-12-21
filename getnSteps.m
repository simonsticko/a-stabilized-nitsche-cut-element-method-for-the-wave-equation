%Returns the timestep and the number of steps that should be taken
%according to cfl condition. cfl is the cfl number,c is wave speed, T is
%end time, n is number of gridpoints over the length domainSize.
function[dt,nSteps]=getnSteps(cfl,c,T,n,domainSize)
%approximative mesh size.
h=domainSize./n;
%maximum timestep
dtMax=cfl*h/c/sqrt(2);
%choose nSteps so that dt is smaller than dtMax.
nSteps=ceil(T./dtMax)+1;
dt=T./(nSteps-1);
end