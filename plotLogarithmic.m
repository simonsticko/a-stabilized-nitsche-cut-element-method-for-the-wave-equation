%Plots vector error against h in log-log and plots a linear least
%square fit of log(h) against log(xAxis). In order to determine the
%dependence error=C*h^p.
function[]=plotLogarithmic(h,error,xLabel,yLabel,expectedSlope,color)
logh=log(h);
logError=log(error);
%draw line through minimal error point.
[logMinError,index]=min(logError);
poly=[expectedSlope; logMinError-expectedSlope*logh(index)];

%plot should be a straight line
plot(logh,logError,[color 'o'],logh,polyval(poly,logh),[color '--'],'linewidth',2);
xlabel(['log(' xLabel ')'],'interpreter','latex');
ylabel(['log(' yLabel ')'],'interpreter','latex');
end