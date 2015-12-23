function[]=processtimestep(ddate)
close all;
datafile=['results/timestep/' ddate '/savedData'];
load(datafile)
% plotEigs(h,epsilon,eigNoStab,[resultdir, 'noStab']);
% plotEigs(h,epsilon,eigOnlyM,[resultdir, 'onlyM']);
plotEigs(h,epsilon,eigOnlyA,[resultdir, 'onlyA']);
plotEigs(h,epsilon,eigStab,[resultdir, 'stab']);
end

function[]=plotEigs(h,epsilon,eigenvalues,savename)
fig=figure();
linewidth=2;
markersize=8;
hMat=h*ones(size(epsilon));
linespecs={':d','-.s','--^'};
eigScaled=hMat.^2.*eigenvalues;
hold on;
for i=1:size(eigenvalues,1)
    plot(log10(epsilon),log10(eigScaled(i,:)),linespecs{i},'linewidth',linewidth,...
        'markersize',markersize);
end
hold off;
%Fix legend
formatfunction=@(s) ['h=' num2str(s,'%1.1e')];
llegend=arrayfun(formatfunction, h, 'unif', 0);
legend(llegend,'orientation','horizontal','location','NorthOutside');
xlabel('$\log_{10}(\epsilon)$','interpreter','latex');
ylabel(['$\log_{10}(h^2\lambda_{\max})$'],'interpreter','latex');
resizeFig(14);
fixPaperSize();
saveas(fig,[savename ,'.pdf'],'pdf');
saveas(fig,[savename ,'.fig'],'fig');
end
