%Assembles the matrix J (see Burman Hansbo page 330 and notes).
function[J]=assembleOverAssociatedFaces(cutMesh)
%Theese are the faces known as FG in article by burman hanso
faces=cutMesh.assocFaces;
X=cutMesh.dt.Points;
nNodes=size(X,1);
J=sparse(nNodes,nNodes);
nFaces=size(faces,1);
for j=1:nFaces
    face=faces(j,:);
    %The face as a vector.
    dX=X(face(1),:)-X(face(2),:);
    %length of this face.
    L=norm(dX);
    normal=[dX(2) -dX(1)]/L;
    %which triangles in the connectiviy list this face is attached to.
    attachments=cutMesh.dt.edgeAttachments(face);
    %convert from cell to vector.
    attachments=attachments{:};
    %The elements connected to this edge, K and Kprime=Kp:
    K=cutMesh.dt.ConnectivityList(attachments(1),:);
    Kp=cutMesh.dt.ConnectivityList(attachments(2),:);
    [gradJump,loc2glb]=gradJumpTwoTriangles(X(K,:),X(Kp,:),K,Kp,normal);
    %Form Jloc from outer product.
    Jloc=L*(gradJump*gradJump');
    J(loc2glb,loc2glb)=J(loc2glb,loc2glb)+cutMesh.h*Jloc;
end
end

%Caclulates the gradient jump term between two triangles K and Kp.
%normal is the unit normal of this face. Output is the gradient jump
%and the 4-by-4 local to global mapping. Note that the gradient jump 
%is 4-by-1.
function[gradJump,loc2glb]=gradJumpTwoTriangles(XK,XKp,K,Kp,normal)
%Local to global mapping should be 4-by-1, since the two elements have
%4 nodes together. Create a local to global mapping and keep track of
%which elements in K and Kp goes where in this loc2glb.
[loc2glb,~,indloc2glb]=unique([K Kp]);
K2loc=indloc2glb(1:3);
Kp2loc=indloc2glb(4:6);
gradJump=gradientJump(XK,XKp,K2loc,Kp2loc,normal);
end

%Returns the gradient jump between 2 triangles K and Kp (K-prime) the
%inparameters are the coordinates the local to global mapping of the
%triangles and the normal of the face that they share.
function[gradJump]=gradientJump(XK,XKp,K2loc,Kp2loc,normal)
%compute gradients.
[~,b,c]=getPhi(XK(:,1),XK(:,2));
[~,bp,cp]=getPhi(XKp(:,1),XKp(:,2));
%Complete gradients on K and Kp in local coordinates, that is "in loc2glb".
gradOnK=zeros(4,2);
gradOnKp=zeros(4,2);
gradOnK(K2loc,:)=[b c];
gradOnKp(Kp2loc,:)=[bp cp];
gradJump=rowWiseScalarProd(gradOnK,normal)-...
    rowWiseScalarProd(gradOnKp,normal);
end
