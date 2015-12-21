%Resizes and saves figure.
function[]=plotAndSave(h,error,yLab,path,expectedSlope)
fig=figure();
name=inputname(2);
color='b';
plotLogarithmic(h,error,'h',yLab,expectedSlope,color);
resizeFig();
fixPaperSize();
saveas(figure(fig),[path name '.fig'],'fig');
saveas(figure(fig),[path name '.pdf'],'pdf');
end