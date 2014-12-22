%Helper function for assemble, assembles over the outer edges of the
%background mesh. See function assemble for an description of the notation.
function[a,L]=assembleOverOuter(cutMesh,gD,dirichlet,gammaD,gammaN)
edges=cutMesh.dt.freeBoundary;
nNodes=size(cutMesh.dt.Points,1);
C=sparse(nNodes,nNodes);
D=sparse(nNodes,nNodes);
E=sparse(nNodes,nNodes);
G=zeros(nNodes,1);
H=zeros(nNodes,1);
for j=1:size(edges,1)
    edge=edges(j,:);
    Xedge=cutMesh.dt.Points(edge,:);
    %Want to now which triangle this edge is attached to.
    attachment=cutMesh.dt.edgeAttachments(edge);
    attachment=attachment{:};
    %Triangle which this cell belongs to.
    loc2glb=cutMesh.dt.ConnectivityList(attachment,:);
    Xelem=cutMesh.dt.Points(loc2glb,:);
    [a,b,c]=getPhi(Xelem(:,1),Xelem(:,2));
    phiAB=a*ones(1,2)+b*Xedge(:,1)'+c*Xedge(:,2)';
    [CLoc,DLoc,ELoc,GLoc,HLoc]=get1DQuantities(Xedge,a,b,c,phiAB,gD);
    %Map local 2 global.
    C(loc2glb,loc2glb)=C(loc2glb,loc2glb)+CLoc;
    D(loc2glb,loc2glb)=D(loc2glb,loc2glb)+DLoc;
    E(loc2glb,loc2glb)=E(loc2glb,loc2glb)+ELoc;
    G(loc2glb)=G(loc2glb)+GLoc;
    H(loc2glb)=H(loc2glb)+HLoc;
end
h=cutMesh.h;
L=(gammaD/h*H-G)*dirichlet+(H+gammaN*h*G)*(~dirichlet);
%This is the terms in a which gives a contribution on the outer boundary.
a=(-C+gammaD/h*D)*dirichlet+(gammaN*h*E)*(~dirichlet);
end