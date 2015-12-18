function[]=processtimestep()
clear;
close all;
load timestep
ploteigs(h,epsilon,stabilized);
end

function[]=ploteigs(h,epsilon,eigenvalues)
linewidth=2;
hMat=h*ones(size(epsilon));
fig=figure();
plot(h,log(hMat.^2.*eigenvalues),'linewidth',linewidth);

end
