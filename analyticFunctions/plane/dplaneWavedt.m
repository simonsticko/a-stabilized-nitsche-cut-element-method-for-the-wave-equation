function[dudt]=dplaneWavedt(x,y,t,k,c)
ssize=size(x);
x=x(:);
y=y(:);
direction=[1,0];
kx=k*(direction*[x';y'])';
omega=k*c;
dudt=real(-1i*k*c*exp(1i*(kx-omega*t)));
dudt=reshape(dudt,ssize);