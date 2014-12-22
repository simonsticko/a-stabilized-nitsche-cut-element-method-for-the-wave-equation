function[r,theta,nTerms,ssize]=softCommon(X,Y)
nTerms=50;%145;
%Save size and rehshape to vectors.
ssize=size(X);
X=X(:);
Y=Y(:);
%Cylinder coordinates.
r = sqrt(Y.^2+X.^2);
theta = atan2(Y,X);
end