%Returns local lagrange multiplier taken over the incoming element.
%this is the local matrix for the matrix with
%lagrange_j=(phi_j,1).
function[lagrangeLoc]=getlagrangeLoc(a,b,c,aarea,Xelem)
%Midpoint on each triangle edge.
Xmids=.5*(Xelem+circshift(Xelem,1));
nVertices=size(Xelem,1);
phiMids=a*ones(1,nVertices)+b*Xmids(:,1)'+c*Xmids(:,2)';
lagrangeLoc=aarea/nVertices*sum(phiMids,2);
end