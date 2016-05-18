%Creates a avi movie file from the saved file with filename loadName,
%the movie is taken from angle view with nFrames number of frames.
function[]=playWave(loadName,vview,nFrames,axisEqual)
if(nargin<2)
    vview=.5*[1 1 1];
    axisEqual=false;
end
videoName=loadName;
[fig1,videoObj]=setupForPlotting(videoName);
load([loadName '.mat']);
%Avoid plotting all calculated values
    if(nargin<3)
        nPlotsPerPeriod=80;
        nFrames=min(nPlotsPerPeriod*nPeriods,length(t));
    end
toPlot=floor(linspace(1,length(t),min(nFrames,length(t))));
% uAnaly=@(t0) besselWave(dTri.Points(relevant,1),dTri.Points(relevant,2),t0,R,c,besselZero);
height=max(max(abs(u-uConst)));
zLim=[uConst-height, uConst+height];
aaxis=[xLim,yLim,zLim];
ccaxis=[uConst-height, uConst+height];
dudtHeight=max(max(abs(dudt)));
dudtaxis=[xLim,yLim,[-dudtHeight,dudtHeight]];
ccaxisdudt=[-dudtHeight,dudtHeight];
set(gcf,'color','w');
for j=1:length(toPlot)
    clf(fig1);
    subplot(1,2,1);
    surfWrap(cutMesh,u(:,toPlot(j)),uConst);
    axis(aaxis);
    if(axisEqual)
        axis equal;
    end
    axis off;
    caxis(ccaxis);
    view(vview);
    title(['t=' sprintf('%.2f',t(toPlot(j)))],'FontSize',52)
    subplot(1,2,2);
    surfWrap(cutMesh,dudt(:,toPlot(j)),0);
    axis(dudtaxis);
    if(axisEqual)
        axis equal;
    end
    axis off;
    caxis(ccaxisdudt);
    view(vview);
    frame=getframe(fig1);
    videoObj.writeVideo(frame);
end
close(videoObj);
end

function[fig1,videoObj]=setupForPlotting(videoName)
% close all;
fig1=figure(1);
a = axes('Parent',fig1);
axis(a,'tight');
set(a,'nextplot','replacechildren');
pos=get(0,'Screensize');
set(fig1, 'Position', [pos(1:2) .5*pos(3:4)]);
videoObj =  VideoWriter([videoName '.avi']);
open(videoObj);
set(fig1,'render','zbuffer');
end