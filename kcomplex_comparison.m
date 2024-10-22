clear all
close all
clc


set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultLineLineWidth',1.5)

fsize = 22;

tcl = tiledlayout(1,3);

b = 0.2;
load(join(['/media/samantha/My Passport/these_sillage_plots/R_kcomplex/', num2str(round(b, 2)), '.mat']));

ks_real = linspace(0.001,pi,100);
ks_imag = linspace(-0.5*pi, 0.5*pi, 100);
ii=1;
nexttile()
pcolor(ks_real, ks_imag, 20*log(abs(reflection)).')

text(0.1, 1.35, '(a)')

hold off
cmap = cmocean('balance');
set(gca, 'colormap', cmap);
shading interp
set(gca, 'clim', [-50 50]);
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})

set(gca,'YTick',-pi/4:pi/4:pi/4) 
set(gca, 'Ylim', [-pi/4 pi/4]);
if ii==1
	set(gca,'YTickLabel',{'$-\pi/4$', '$0$', '$\pi/4$'})
	%set(gca,'YTickLabel',{'$-\pi/2$','$-\pi/4$', '$0$', '$\pi/4$','$\pi/2$'})
else
	set(gca,'YTickLabel',{'','', ''})
	%set(gca,'YTickLabel',{'','', '', '',''})
end
xlabel('Re($ka$)')
if ii==1, ylabel('Im($ka$)'), end
%title(join(['$b/a$ = ', num2str(round(b, 2))]))
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
title('Multimodal method',  'FontSize', 24)

ii=2;
nexttile()

%% specular reflection

theta = 35*pi/180;
L = 1;
for bb=1:length(ks_real)
	for dd=1:length(ks_imag)
		k = ks_real(bb) + 1j*ks_imag(dd);
		%kx = k*sin(theta); % kx0
		ky = k*cos(theta); %ky0
		
		num = 1j.*ky - k.*b.*tan(k.*L);
		den = 1j.*ky + k.*b.*tan(k.*L);
		Rspecular(bb, dd) = num./den;
	end
end


pcolor(ks_real, ks_imag, 20*log(abs(Rspecular)).')

text(0.1, 1.35, '(b)')
hold off
cmap = cmocean('balance');
set(gca, 'colormap', cmap);
shading interp
set(gca, 'clim', [-50 50]);
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})
set(gca,'YTick',-pi/4:pi/4:pi/4) 
set(gca, 'Ylim', [-pi/4 pi/4]);
if ii==1
	set(gca,'YTickLabel',{'$-\pi/4$', '$0$', '$\pi/4$'})
	%set(gca,'YTickLabel',{'$-\pi/2$','$-\pi/4$', '$0$', '$\pi/4$','$\pi/2$'})
else
	set(gca,'YTickLabel',{'','', ''})
	%set(gca,'YTickLabel',{'','', '', '',''})
end
xlabel('Re($ka$)')
if ii==1, ylabel('Im($ka$)'), end
%title(join(['$b/a$ = ', num2str(round(b, 2))]))
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
title('Effective boundary condition',  'FontSize', 24)

ii=3;
nexttile()

%% specular reflection


Rdots = abs(20*log(abs(Rspecular)))>60;

%pcolor(ks_real, ks_imag, 20*log(abs(Rdots)).')
contour(ks_real, ks_imag,  20*log(abs(Rspecular)).', 20)

hold on
plot(pi/2, mean([ks_imag(58), ks_imag(59)]), 'ko')
hold on
plot(pi/2, mean([ks_imag(42), ks_imag(43)]), 'ko')
plot(ks_real, ones([1,length(ks_real)])*0, 'k--')

x1p = pi/2;
y1p = mean([ks_imag(58), ks_imag(59)]);
x2p = x1p;
y2p = 0;
arh2 = annotation('arrow', 'HeadStyle', 'plain', 'LineWidth', 1,'HeadLength', 5, 'HeadWidth', 5, 'color', 'k');
arh2.Units = 'normalized';
ax = gca;
arh2.Parent = ax;
arh2.Position = [x1p, y1p, x2p-x1p, y2p-y1p];
%arh2.Position = [0.3, 0.3, sin(theta), cos(theta)];
arh2.Color = 'k';

x1p = pi/2;
y1p = mean([ks_imag(42), ks_imag(43)]);
x2p = x1p;
y2p = y1p-mean([ks_imag(58), ks_imag(59)]);
arh2 = annotation('arrow', 'HeadStyle', 'plain', 'LineWidth', 1,'HeadLength', 5, 'HeadWidth', 5, 'color', 'k');
arh2.Units = 'normalized';
ax = gca;
arh2.Parent = ax;
arh2.Position = [x1p, y1p, x2p-x1p, y2p-y1p];
%arh2.Position = [0.3, 0.3, sin(theta), cos(theta)];
arh2.Color = 'k';

text(pi/2 + 0.15, mean([ks_imag(58), ks_imag(59)]), 'zero')
text(pi/2 + 0.15, mean([ks_imag(42), ks_imag(43)]), 'pole')

text(0.1, 1.35, '(c)')


hold off
cmap = cmocean('balance');
set(gca, 'colormap', cmap);
shading interp
set(gca, 'clim', [-50 50]);
%set(gca, 'Ylim', [-pi/2 pi/2]);
set(gca, 'Xlim', [0 pi]);
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})
set(gca,'YTick',-pi/4:pi/4:pi/4) 
set(gca, 'Ylim', [-pi/4 pi/4]);
if ii==1
	set(gca,'YTickLabel',{'$-\pi/4$', '$0$', '$\pi/4$'})
	%set(gca,'YTickLabel',{'$-\pi/2$','$-\pi/4$', '$0$', '$\pi/4$','$\pi/2$'})
else
	set(gca,'YTickLabel',{'','', ''})
	%set(gca,'YTickLabel',{'','', '', '',''})
end
xlabel('Re($ka$)')
if ii==1, ylabel('Im($ka$)'), end
title('Lossy system', 'FontSize', 24)
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)


cb = colorbar('TickLabelInterpreter', 'latex');
cb.Layout.Tile = 'east'; % Assign colorbar location

saveas(tcl,'/home/samantha/Dropbox/PhD/these/figures/wake/kcomplex_absorption.png')
