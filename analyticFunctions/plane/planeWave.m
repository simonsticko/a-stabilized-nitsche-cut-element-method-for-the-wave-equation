function[u]=planeWave(x,y,t,k,c)
ssize=size(x);
x=x(:);
y=y(:);
direction=[1,0];
kx=k*(direction*[x';y'])';
omega=k*c;
u=real(exp(1i*(kx-omega*t)));
u=reshape(u,ssize);