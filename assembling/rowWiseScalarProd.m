%Row wise scalar product, if Z is a n-by-2 matrix and X is 2-by-1 sprod is
%the n-by-1 vector where sprod_j=dot(Z(j,:),X)
function[sprod]=rowWiseScalarProd(Z,X)
sprod=sum(Z.*kron(ones(size(Z,1),1),X),2);
end