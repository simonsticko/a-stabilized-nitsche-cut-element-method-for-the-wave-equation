function[graddudt]=dgradPlaneWavedt(x,y,t,k,c)
x=x(:);
y=y(:);
direction=[1,0];
kx=k*(direction*[x';y'])';
omega=k*c;
graddudt=real((-1i*k*c)*(1i*k)*exp(1i*(kx-omega*t)))*direction;