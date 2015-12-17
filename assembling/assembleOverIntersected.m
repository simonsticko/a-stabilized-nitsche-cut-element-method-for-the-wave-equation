%Assembles over the triangles that are intersected by the boundary.
%helper function for assemble. see text in assemble for description. Of
%some of the parameters.
%G and H are defines as follows:
%Gj=<g,phi_j>
%Hj=<g,grad_phi_j>
%and are used for the boundary enforcing terms in L:
%if using dirichlet boundary conditions we have
%L=bLoad+(gammaD/h*H-G);
%and if using neumann bc we have
%L=bLoad+(H+gammaN*h*G);
function[m,a,L,lagrange]=...
    assembleOverIntersected(cutMesh,f,g,dirichlet,gammaD,gammaN)
X=cutMesh.dt.Points(:,:);
nNodes=size(X,1);
%Preallocate.
m=sparse(nNodes,nNodes);
lagrange=zeros(nNodes,1);
B=sparse(nNodes,nNodes);
C=sparse(nNodes,nNodes);
D=sparse(nNodes,nNodes);
E=sparse(nNodes,nNodes);
G=zeros(nNodes,1);
H=zeros(nNodes,1);
bLoad=zeros(nNodes,1);
%Loop over all intersected elements
tri=cutMesh.getTriIntersected();
for j=1:size(tri,1)
    loc2glb=tri(j,:);
    Xelem=X(loc2glb,:);
    %Intersected element, note that this can have 3 or for "vertices".
    Xcut=cutMesh.Xcut{j};
    [a,b,c]=getPhi(Xelem(:,1),Xelem(:,2));
    aarea=polyarea(Xcut(:,1),Xcut(:,2));
    nVertices=size(Xcut,1);
    %Evaluate phij at the ndoes of the cut element each row is phi evaluated
    %at the vertices of the cut element.
    phiN=a*ones(1,nVertices)+b*Xcut(:,1)'+c*Xcut(:,2)';
    fN=ones(3,1)*f(Xcut(:,1),Xcut(:,2))';
    %Quadrature see eq 3.57 in Larson Bengzon.
    bLoadLoc=sum(fN.*phiN,2)*aarea/nVertices;
    BLoc=aarea*(b*b'+c*c');
    [mLoc]=getmLoc(phiN,nVertices,aarea);
    [lagrangeLoc]=getlagrangeLoc(a,b,c,aarea,Xcut);
    [CLoc,DLoc,ELoc,GLoc,HLoc]=get1DQuantities(Xcut,...
        a,b,c,phiN(:,1:2),g);
    %Map local to global:
    m(loc2glb,loc2glb)=m(loc2glb,loc2glb)+mLoc;
    bLoad(loc2glb)=bLoad(loc2glb)+bLoadLoc;
    lagrange(loc2glb)=lagrange(loc2glb)+lagrangeLoc;
    B(loc2glb,loc2glb)=B(loc2glb,loc2glb)+BLoc;
    C(loc2glb,loc2glb)=C(loc2glb,loc2glb)+CLoc;
    D(loc2glb,loc2glb)=D(loc2glb,loc2glb)+DLoc;
    E(loc2glb,loc2glb)=E(loc2glb,loc2glb)+ELoc;
    G(loc2glb)=G(loc2glb)+GLoc;
    H(loc2glb)=H(loc2glb)+HLoc;
end
h=cutMesh.h;
L=bLoad+(gammaD/h*H-G)*dirichlet+(H+gammaN*h*G)*(~dirichlet);
a=B+(-C+gammaD/h*D)*dirichlet+(gammaN*h*E)*(~dirichlet);
end

