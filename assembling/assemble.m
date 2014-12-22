%Assembles the system matrices.
%M- mass matrix with stabilization.
%A- RHS matrix from burman-hansbo.
%lagrange is a lagrange multiplier which can be used for a pure neumann
%problem. lagrange_j=(phi_j,1). This has not been tested though.
%We use the following notation for the sub-matrices of A:
%Bij=(phi_i,phi_j)
%Cij=<phi_i,n*grad_phi_j>+<n*grad_phi_i,phi_j>
%Dij=<phi_i,phi_j>
%Eij=<n*grad_phi_i,n*grad_phi_j>
%Where (,) denotes L2-scalar product over the 2d domain and <,> denotes L2
%scalar prodduct over the boundary.
%Furhermore denote bLoad_j=(f,phi_j).
%
%Ending Inside in the name of a matrix means that this is the part of the matrix that
%is assembled over the inside of the mesh. 
%Ending Inter means that this is the part of the matrix that is assembled
%over the elements cut by the boundary.
%Ending out means the part of the matrix that is assembled over the outer
%edges of the background mesh. This part is only relevant when solving an
%outer problem.
%
%a and m denotes mass and stiffness matrix without stabilization.
%
%Inparamters
%cutMesh is an object of type CutMesh
%f is the forcing function of the system
%gcSq is the function in the boundary condition scaled with the 
%wave speed squared: gcSq=c^2*g. Where u=g dirichletInner=true
%and n*grad_u=g if dirichletInnter=false. 
%gOutcSq and dirichletOuter works in the same way and is only needed when
%solving and outer problem, for the outer problem g will be the inner
%boundary condition, for an inner problem this is the only boundary
%condition.
function[M,A,L,LOut,lagrange]=assemble(cutMesh,f,gcSq,dirichletInner,gOutcSq,dirichletOuter)
%Penallty constants.
gamma1=.5;
gammaM=.25;
gammaD=5;
gammaN=5;
if(nargin<5)
   gOutcSq=@(x,y) zeros(size(x));
   dirichletOuter=true; 
end
%small m denotes mass matrix without stabilization.
[mInside,aInside,LInside,lagrangeInside]=assembleOverInside(cutMesh,f);
[mInter,aInter,LInter,lagrangeInter]=...
    assembleOverIntersected(cutMesh,f,gcSq,dirichletInner,gammaD,gammaN);
%Stabilizing term.
J=assembleOverAssociatedFaces(cutMesh);
%Boundary condition for outer problem:
[aOut,LOut]=assembleOverOuter(cutMesh,gOutcSq,dirichletOuter,gammaD,gammaN);
%get nodes which is relevant for the calculation.
relevant=cutMesh.relevant;
%Add the contributions from the inner, intersected and outer assembling.
m=mInside(relevant,relevant)+mInter(relevant,relevant);
lagrange=lagrangeInside(relevant)+lagrangeInter(relevant);
a=aInter(relevant,relevant)+aInside(relevant,relevant)+aOut(relevant,relevant);
L=LInside(relevant)+LInter(relevant);
%Loutneeds to be a separate output since the outer problem has also a
%temporal dependence which one can separate multiplicatively.
LOut=LOut(relevant);
J=J(relevant,relevant);
%add stabilization.
h=cutMesh.h;
M=m+gammaM*h^2*J;
A=a+gamma1*J;
end

%Assembles over the inside of the domain.
function[m,B,bLoad,lagrange]=assembleOverInside(cutMesh,f)
X=cutMesh.dt.Points();
%Total numer of nodes in the system.
nNodes=size(cutMesh.dt.Points(),1);
%Preallocate.
m=sparse(nNodes,nNodes);
B=sparse(nNodes,nNodes);
bLoad=zeros(nNodes,1);
lagrange=zeros(nNodes,1);
%Each inside triangle has 3 vertices, since they are not cut.
nVertices=3;
%Loop over all intersected elements
tri=cutMesh.getTriInside();
for j=1:size(tri,1)
    loc2glb=tri(j,:);
    Xelem=X(loc2glb,:);
    %Two a in aarea because area is a function;
    [a,b,c,aarea]=getPhi(Xelem(:,1),Xelem(:,2));
    BLoc=aarea*(b*b'+c*c');
    %This quadrature comes from eq 3.57 in Larson Bengzon by noting that
    %phi_i is zero on the other nodes of the element.
    bLoadLoc=f(Xelem(:,1),Xelem(:,2))*aarea/3;
    phiN=a*ones(1,nVertices)+b*Xelem(:,1)'+c*Xelem(:,2)';
    MLoc=getmLoc(phiN,nVertices,aarea);
    [lagrangeLoc]=getlagrangeLoc(a,b,c,aarea,Xelem);
    %Map local to global:
    B(loc2glb,loc2glb)=B(loc2glb,loc2glb)+BLoc;
    bLoad(loc2glb)=bLoad(loc2glb)+bLoadLoc;
    m(loc2glb,loc2glb)=m(loc2glb,loc2glb)+MLoc;
    lagrange(loc2glb)=lagrange(loc2glb)+lagrangeLoc;
end
end
