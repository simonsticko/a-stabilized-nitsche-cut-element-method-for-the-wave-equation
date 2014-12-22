%Returns the gradient of the n'th "mode" solution to the circular membrane
%problem.
function[gradu,T]=gradBesselWave(x,y,t,R,c,n)
r=sqrt(x.^2+y.^2);
rHat=[x./r,y./r];
%zeroth order bessel function zeros.
z=[2.4048 5.5201 8.6537 11.7915 14.9309];
zn=z(n);
% z0=fzero(j0,2.4);
T=2*pi*R/zn/c;
absgradu=-zn/R*besselj(1,r/R*zn).*(r<=R)*cos(2*pi*t/T);
gradu=[rHat(:,1).*absgradu, rHat(:,2).*absgradu];
end