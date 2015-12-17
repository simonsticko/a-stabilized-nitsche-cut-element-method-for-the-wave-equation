clear;
n=100;
nT=1E4;
x=linspace(-2.5,2.5,n);
y=linspace(-2.5,2.5,n);
[X,Y]=meshgrid(x,y);
out=sqrt(X.^2+Y.^2)>1;
gamma=pi;
F=zeros(size(out));
F(out)=softCylSpat(X(out),Y(out),gamma);
Xu=X(:);
Yu=Y(:);
outu=out(:);
gradF=zeros(length(outu),2);
gradF(outu,:)=gradSoftCylSpat(Xu(outu),Yu(outu),gamma);
u=@(t) real(F*exp(1i*gamma*t));
gradu=@(t) real(gradF*exp(1i*gamma*t));
t=linspace(0,100,nT)';
for j=1:length(t)
    subplot(1,2,1);
    surf(X,Y,u(t(j)));
    view(2);
    subplot(1,2,2);
    gradut=gradu(t(j));
    quiver(X(:),Y(:),gradut(:,1),gradut(:,2));
    axis equal;
    pause(0.01);
end