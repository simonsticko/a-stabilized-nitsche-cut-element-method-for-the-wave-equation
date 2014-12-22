%Returns the local mass matrix. phiN is a 3-by-nVertices matrix. Each
%element phi_ij is hat function phi_i evaluated at the vertex j.
%aarea is the area of the element.
function[mLoc]=getmLoc(phiN,nVertices,aarea)
    %Local mass matrix
    mLoc=zeros(3,3);
    for j=1:nVertices
        phij=phiN(:,j);
        mLoc=mLoc+aarea/nVertices*(phij*phij');
    end
end