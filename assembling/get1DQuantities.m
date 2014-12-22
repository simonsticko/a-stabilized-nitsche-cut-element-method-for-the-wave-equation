%This function is helper for assembleOverIntersected.
%Returns the one dimensional quantities integrated over Xcut(1:2,:) which
%is the edge part of the intersected triangle.
%This edge will be the straight line integrated from XA to XB were
%XA=Xcut(1,:) and XB=Xcut(2,:);
function[CLoc,DLoc,ELoc,GLoc,HLoc]=get1DQuantities(Xcut,a,b,c,phiAB,g)
%get normal to this edge, L is the length of XA-XB.
[normal,L,XA,XB]=getNormal(Xcut);
gradPhi=[b c];
gradPhiDotn=rowWiseScalarProd(gradPhi,normal);
[CLoc,DLoc,ELoc]=getCLocDLocELoc(gradPhi,gradPhiDotn,a,XA,...
    XB,L);
%Boundary enforcing terms see definition in assembleOverIntersected.
[GLoc,HLoc]=getGandH(g,XA,XB,phiAB,gradPhiDotn,L);
end

%See assembleOverIntersected for definition of G and H.
function[G,H]=getGandH(g,XA,XB,phiAB,gradPhiDotn,L)
gA=g(XA(1),XA(2));
gB=g(XB(1),XB(2));
H=L/2*(gA*phiAB(:,1)+gB*phiAB(:,2));
G=L/2*(gA+gB)*gradPhiDotn;
end

%Returns the local matrices in the stiffness matrix.
function[CLoc,DLoc,ELoc]=getCLocDLocELoc(gradPhi,gradPhiDotn,a,XA,XB,L)
XBA=XB-XA;
%Note that these are not the bc and forcing function.
g=a+rowWiseScalarProd(gradPhi,XA);
f=rowWiseScalarProd(gradPhi,XBA);
%A for assymmetric. Note that this is the outer product.
CLocA=L*gradPhiDotn*(g+.5*f)';
CLoc=CLocA+CLocA';
fg=f*g';
DLoc=L*(g*g'+(fg+fg')/2+(f*f')/3);
ELoc=L*(gradPhiDotn*gradPhiDotn');
end