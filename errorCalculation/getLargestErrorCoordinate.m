%Returns the coordinates of the node in uEnd having largest error. dTri is
%the delauneytriangulation, and uAnalyEnd is the analytic function. XB is
%the boundary polygon and relevant is the logical array defining which
%nodes are relevant for the calculation.
function[XMaxError]=getLargestErrorCoordinate(uEnd,dTri,uAnalyEnd,XB,relevant)
uFullEnd=zeros(size(dTri.Points,1),1);
uFullEnd(relevant)=uEnd;
endErrorFull=abs(uFullEnd-uAnalyEnd(dTri.Points(:,1),dTri.Points(:,2)));
%Set error outside domain to zero.
nodesInside=inpolygon(dTri.Points(:,1),dTri.Points(:,2),XB(:,1),XB(:,2));
endErrorFull(~nodesInside)=0;
[~,indexMax]=max(endErrorFull);
XMaxError=dTri.Points(indexMax,:);
end