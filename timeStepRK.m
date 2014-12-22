%Timesteps the system obtained from discretizing the wave equation, with a
%classical 4th order runge kutta.
function[u,dudt,t]=timeStepRK(u0,du0dt,A,M,L,nSteps,dt,c)
%t must be a row vector
t=dt*(0:(nSteps-1));
%A has dimenstion length^-2 this makes D2 be dimensionless
D2=-c^2*A;
%LU-factorize M for use later
[Llu,Ulu]=lu(M);
% %Solution matrix, intiutivt: v=[u dudt]
v=zeros(2*length(u0),nSteps);
% %inital condition derivative is zero
% v(1:(end/2),1)=u0;
v(:,1)=[u0;du0dt];
dvdt=@(t,w) dvdtD2LU(t,w,D2,Llu,Ulu,L);
for j=1:(nSteps-1)
    tj=t(j);
    vj=v(:,j);
    K1=dt*dvdt(tj,vj);
    K2=dt*dvdt(tj+.5*dt,vj+K1/2);
    K3=dt*dvdt(tj+.5*dt,vj+K2/2);
    K4=dt*dvdt(tj+dt,vj+K3);
    v(:,j+1)=vj+1/6*(K1+2*K2+2*K3+K4);
end
u=v(1:(end/2),:);
dudt=v((end/2)+1:end,:);
end

%Function returning the derivative of our system.
function[dvdtn]=dvdtD2LU(t,vn,D2,Llu,Ulu,L)
nNodes=length(vn)/2;
dvdtn=zeros(size(vn));
dvdtn(1:nNodes)=vn((nNodes+1):end);
dvdtn((nNodes+1):end)=Ulu\(Llu\(D2*vn(1:nNodes)+L(t)));
end
