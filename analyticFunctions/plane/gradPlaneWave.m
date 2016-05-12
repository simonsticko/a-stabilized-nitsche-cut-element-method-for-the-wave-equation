function[gradu]=gradPlaneWave(x,y,t,k,c)
x=x(:);
y=y(:);
direction=[1,0];
kx=k*(direction*[x';y'])';
omega=k*c;
gradu=real(1i*k*exp(1i*(kx-omega*t)))*direction;