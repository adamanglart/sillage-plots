clear all
close all
clc


set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultLineLineWidth',1.5)

fsize = 20;

bs = 0.1:0.1:0.3;
bs = [0.2];
%tcl = tiledlayout(1,3);
tcl = tiledlayout(1,1);
for ii=1:length(bs)

b = bs(ii);
%b = 0.1;

load(join(['/media/samantha/My Passport/these_sillage_plots/R_kcomplex/', num2str(round(b, 2)), '.mat']));

ks_real = linspace(0.001,pi,100);
ks_imag = linspace(-0.5*pi, 0.5*pi, 100);
nexttile()
%figure;
pcolor(ks_real, ks_imag, 20*log(abs(reflection)).')
%hold on
%contourf(ks_real, ks_imag, 20*log(abs(reflection)).', 'levels', 4)
%plot(ks_real, 0*ks_real, 'k--')
hold off
%dx = 0.2
%xlim([1.46-dx, 1.46+dx])
%ylim([-1, 1])
cmap = cmocean('balance');
set(gca, 'colormap', cmap);
shading interp
set(gca, 'clim', [-50 50]);
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})
set(gca,'YTick',-pi/2:pi/4:pi/2) 
if ii==1
	set(gca,'YTickLabel',{'$-\pi/2$','$-\pi/4$', '$0$', '$\pi/4$','$\pi/2$'})
else
	set(gca,'YTickLabel',{'','', '', '',''})
end
xlabel('Re($ka$)')
if ii==1, ylabel('Im($ka$)'), end
title(join(['$b/a$ = ', num2str(round(b, 2))]))
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)

end

cb = colorbar('TickLabelInterpreter', 'latex');
cb.Layout.Tile = 'east'; % Assign colorbar location
saveas(tcl,'/home/samantha/Dropbox/PhD/these/figures/wake/kcomplex_absorption.png')
