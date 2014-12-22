%Returns the RHS that should be used to get the ritz projection. That is
%Au is the vector Au_j=a(u,phi_j). u0Analy and gradu0Analy is the analytic
%function and gradient of this analytic function.
%
%Matrix endings Inside, Inter and Out is described in assemble.
%
%Endings D and B on Au denotes the terms of Au taken over the boundary and
%the domain respectively. That is:
%AuD_j=(grad_u,grad_phi_j) and AuB denotes the remaining terms, see
%getAuBloc for a description.
function[Au]=getAuRitz(cutMesh,u0Analy,gradu0Analy,dirichletInner,dirichletOuter)
if(nargin<5)
    dirichletOuter=true;
end
%Penalty Constants.
gammaD=5;
gammaN=5;
AuInside=getAuInside(cutMesh,gradu0Analy);
[AuInter]=getAuInters(cutMesh,gradu0Analy,u0Analy,gammaD,gammaN,dirichletInner);
[AuOut]=getAuBOuter(cutMesh,gradu0Analy,u0Analy,gammaD,gammaN,dirichletOuter);
Au=AuInside+AuInter+AuOut;
Au=Au(cutMesh.relevant);
end

%Returns contribution of Au over the inside triangles.
function[AuD]=getAuInside(cutMesh,gradu0Analy)
X=cutMesh.dt.Points;
nNodes=size(X,1);
AuD=zeros(nNodes,1);
tri=cutMesh.getTriInside();
for j=1:size(tri,1)
    loc2glb=tri(j,:);
    Xelem=X(loc2glb,:);
    [~,b,c,aarea]=getPhi(Xelem(:,1),Xelem(:,2));
    AuDLoc=getAuDLoc(Xelem,aarea,b,c,gradu0Analy);
    AuD(loc2glb)=AuD(loc2glb)+AuDLoc;
end
end

%Assembles Au over the outer edge. Needed for the outer problem.
function[Au]=getAuBOuter(cutMesh,gradu0Analy,u0Analy,gammaD,gammaN,dirichlet)
edges=cutMesh.dt.freeBoundary;
nNodes=size(cutMesh.dt.Points,1);
Au=zeros(nNodes,1);
for j=1:size(edges,1)
    edge=edges(j,:);
    Xedge=cutMesh.dt.Points(edge,:);
    attachment=cutMesh.dt.edgeAttachments(edge);
    attachment=attachment{:};
    %Triangle which this cell belongs to.
    loc2glb=cutMesh.dt.ConnectivityList(attachment,:);
    Xelem=cutMesh.dt.Points(loc2glb,:);
    [a,b,c]=getPhi(Xelem(:,1),Xelem(:,2));
    AuBLoc=getAuBLoc(u0Analy,gradu0Analy,Xedge,a,b,c,gammaD,gammaN,...
        cutMesh.h,dirichlet);
    %Map local 2 global.
    Au(loc2glb)=Au(loc2glb)+AuBLoc;
end
end

%Assembled AuD and AuB over the
function[Au]=getAuInters(cutMesh,gradu0Analy,u0Analy,gammaD,gammaN,dirichlet)
X=cutMesh.dt.Points;
nNodes=size(X,1);
Au=zeros(nNodes,1);
tri=cutMesh.getTriIntersected();
for j=1:size(tri,1)
    loc2glb=tri(j,:);
    Xelem=X(loc2glb,:);
    Xcut=cutMesh.Xcut{j};
    [a,b,c]=getPhi(Xelem(:,1),Xelem(:,2));
    aarea=polyarea(Xcut(:,1),Xcut(:,2));
    %Calculate AuD.
    AuDLoc=getAuDLoc(Xcut,aarea,b,c,gradu0Analy);
    AuBLoc=getAuBLoc(u0Analy,gradu0Analy,Xcut,a,b,c,gammaD,gammaN,...
        cutMesh.h,dirichlet);
    %local 2 global.
    Au(loc2glb)=Au(loc2glb)+AuDLoc+AuBLoc;
end
end

%Returns the local AuD term integrated over Xintegration.
%XIntegration is the points that define the domain we integrate over, that
%is Xcut or Xelem.
function[AuDLoc]=getAuDLoc(XIntegration,aarea,b,c,gradu0Analy)
gradU=gradu0Analy(XIntegration(:,1),XIntegration(:,2));
intgradU=aarea/size(XIntegration,1)*sum(gradU,1);
AuDLoc=b*intgradU(1)+c*intgradU(2);
end

%Returns the boundary terms AuB.
function[AuBLoc]=...
    getAuBLoc(u0Analy,gradu0Analy,Xcut,a,b,c,gammaD,gammaN,h,dirichlet)
%Evaluate phij at the ndoes of the cut element each row is phi evaluated
%at the edges of the cut element.
gradPhi=[b c];
%get the outer normal to the boundary.
[normal,L,XA,XB]=getNormal(Xcut);
phiA=a+b*XA(1)+c*XA(2);
phiB=a+b*XB(1)+c*XB(2);
gradPhiDotn=rowWiseScalarProd(gradPhi,normal);
nGradU=@(X) dot(gradu0Analy(X(1),X(2)),normal);
uAn=@(X) u0Analy(X(1),X(2));
%We use the following notation:
%au2_j=<n*grad_u,phi_j>
%au3_j=<u,grad_phi_j>
%au4_j=<u,phi_j>
%au5_j=<n*grad_u,n*grad_phi>
%so that when using dirichlet bc we have:
%AuB=-au2-au3+gammaD/h*au4
%and when using neumman bc we have
%AuB=gammaN*h*au5
au2=L/2*(nGradU(XA)*phiA+nGradU(XB)*phiB);
au3=L/2*(uAn(XA)+uAn(XB))*gradPhiDotn;
au4=L/2*(uAn(XA)*phiA+uAn(XB)*phiB);
au5=L/2*((nGradU(XA)+nGradU(XB))*gradPhiDotn);
AuBLoc=dirichlet*(-au2-au3+gammaD/h*au4)+(~dirichlet)*gammaN*h*au5;
end
