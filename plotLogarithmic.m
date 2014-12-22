%Plots vector error against h in log-log and plots a linear least
%square fit of log(h) against log(xAxis). In order to determine the
%dependence error=C*h^p.
function[p]=plotLogarithmic(h,error,xLabel,yLabel,dydxString,color)
logh=log(h);
logError=log(error);
poly=polyfit(logh,logError,1);
p=poly(1);
%plot should be a straight line
plot(logh,logError,[color 'o'],logh,polyval(poly,logh),[color '--'],'linewidth',2);
xlabel(['log(' xLabel ')'],'interpreter','latex');
ylabel(['log(' yLabel ')'],'interpreter','latex');
%write polynomial in figure, determine where to put the text.
xText=min(logh)+.1*(max(logh)-min(logh));
yText=min(logError)+.9*(max(logError)-min(logError));
text(xText,yText,[dydxString '=' sprintf('%.2f',poly(1))],...
    'fontsize',18,'interpreter','latex');
end