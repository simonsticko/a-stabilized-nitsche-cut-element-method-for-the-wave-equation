%Returns the the n'th "mode" solution to the circular membrane
%problem and the period of the oscillations.
function[u,T]=besselWave(x,y,t,R,c,n,uConst)
r=sqrt(x.^2+y.^2);
%zeroth order bessel function zeros.
z=[2.4048 5.5201 8.6537 11.7915 14.9309];
zn=z(n);
j0=@(z) besselj(0,z);
T=2*pi*R/zn/c;
u=uConst+(j0(r/R*zn).*(r<=R))*cos(2*pi*t/T);
end