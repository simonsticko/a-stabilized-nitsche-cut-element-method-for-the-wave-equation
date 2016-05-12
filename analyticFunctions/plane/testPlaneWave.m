clear;
close all;
n=50;
x=linspace(-1.5,1.5,n);
[X,Y]=meshgrid(x,x);
k=3;
c=1;
time=0;
figure;
surf(X,Y,planeWave(X,Y,time,k,c));
figure;
surf(X,Y,dplaneWavedt(X,Y,time,k,c))
figure;
grad=gradPlaneWave(X(:),Y(:),time,k,c);
quiver(X(:),Y(:),grad(:,1),grad(:,2))
figure;
graddudt=dgradPlaneWavedt(X(:),Y(:),time,k,c);
quiver(X(:),Y(:),graddudt(:,1),graddudt(:,2))