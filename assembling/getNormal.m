%For a given Xcut identifies the outer edges as the first two points of the
%polygon Xcut and returns the normal length of edge and start and end point
%of the edge.
function[normal,L,XA,XB]=getNormal(Xcut)
XA=Xcut(1,:);
XB=Xcut(2,:);
XBA=XB-XA;
%Length of the boundary that cuts this element.
L=norm(XBA);
%unit normal to the boundary segment that cuts this element.
normal=[XBA(2) -XBA(1)]/L;
end