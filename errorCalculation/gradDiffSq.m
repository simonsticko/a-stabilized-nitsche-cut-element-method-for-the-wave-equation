function[diffSq]=gradDiffSq(x,y,gradUAnaly,graduNum)
sizeX=size(x);
x=x(:);
y=y(:);
diffGrad=gradUAnaly(x,y)-graduNum(x,y);
%The 2 means take column wise scalar product.
diffSq=dot(diffGrad,diffGrad,2);
diffSq=reshape(diffSq,sizeX);
end