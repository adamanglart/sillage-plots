clear all
close all
clc

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',1.5)

fsize = 24;

Lx = 10; % length x
Ly = 10; % length y
filling_ratio = 0.5;

%Ly = Ly-1;

R01 = [3,4,-Lx,Lx+0.5,Lx+0.5,-Lx,-Ly,-Ly,0,0]';

n=2*Lx;    % number of plates
%n = 0

P=zeros(n,10);

ns1=blanks(n+1);
ns2=blanks(n+1);
ns3=blanks(n+1);
sf=strings;
ns1(1)='R';
ns2(1)='0';
ns3(1)='1';

h0 = 0.5;
x0 = -Lx -0.5;
%dy = 3*h0;
%dy = 4;
dy = 1;

for jj=1:n
     P(jj,:) = [3,4, x0 + jj, x0 + jj+filling_ratio, x0 + jj+filling_ratio, x0 + jj, -dy, -dy, 0, 0];
    ns1(jj+1)='P';
    temp2=sprintf('%d', fix(jj/10));
    temp3=sprintf('%d', rem(jj,10));
    ns2(jj+1)=temp2;
    ns3(jj+1)=temp3;
    sf=[sf, '-P', temp2 , temp3];
end

P=P';
P=[P;zeros(length(R01)-size(P,1),n)];

geom = [R01,P];
ns=[ns1;ns2;ns3];
sf=['R01',sf];
sf=strjoin(sf,'');
gd = decsg(geom,sf,ns);

model = createpde;
geometryFromEdges(model,gd);

figure(1); 
pdegplot(model,'EdgeLabels','on', 'FaceLabels', 'on'); 
%pdegplot(model,'FaceLabels', 'on'); 
axis equal



bWall = applyBoundaryCondition(model,'neumann','Edge',[1:model.Geometry.NumEdges],'g',0,'q',0);
%bWall = applyBoundaryCondition(model,'dirichlet','Edge',[1:model.Geometry.NumEdges]);

specifyCoefficients(model,'m',0,'d',1,'c',1,'a',0,'f',0, 'Face', 1);

%r = [0,20*pi];
%r = [0,5];
%r = [1.44,1.54];
r = [0, 2];

generateMesh(model,'Hmax',0.2);

figure;
pdemesh(model)

results = solvepdeeig(model,r);

l = results.Eigenvalues;

u = results.Eigenvectors;


N = length(l);


for ii=1:N
	fig = figure(ii);
	pdeplot(model,'XYData',u(:,ii));
	minu = min(real(u(:,ii)));
	maxu = max(real(u(:,ii)));
	
	lim  = max(abs(minu), abs(maxu));


	title(join(['mode ', num2str(ii)]))
	subtitle(join(['$k^2$ = ', num2str(round(l(ii),2))]), 'FontSize',fsize)
	%colormap jet
	cmocean('balance');
	axis equal
	axis off
	ax = gca;
	ax.CLim = [-lim lim];
	%cbar = colorbar;
	%cbar.Position = [0.9, 0.4, 0.0374, 0.23];
	colorbar('off')
	set(fig,'WindowStyle','docked');

end




%% Relacion de dispersion

% medidas del setup (adimensionales) 2022
L = dy;
b = filling_ratio;
a = 1;

% Relacion de dispersion entre q y k
k = linspace(0,pi/(2*L),1000);
q2 = (((b/a)*tan(k*L)).^2 + 1).*(k.^2);
q = sqrt(q2);
%qn = q*b/pi;
qn = q*2*Lx/pi;

% interpolo en los numeros correspondientes a los modos en el grating

% Nmodos = 500;
% qn_interp = 0:1:Nmodos;
% kn_interp = interp1(qn, k, qn_interp);
% 
% figure;
% ax = gca;
% %plot(qn,k, 'Color', '#4169E1', 'DisplayName', 'Dispersion relation');
% plot(qn_interp, kn_interp, '-', 'MarkerFaceColor', 'w', 'MarkerEdgeColor','#4169E1', 'DisplayName', 'Dispersion relation')
% %ylim([min(k) pi/(2*L)])
% xlim([0 Nmodos])
% hold on
% %plot(qn_interp, kn_interp, '-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor','#4169E1')
% 
% % Grafico los modos de superficie que salen de pdetool
% 
% idx_modos_surf = [9, 13, 17, 19, 21, 23, 26, 28, 31, 33, 34, 35, 36, 37, 39, 40, 41];
% k_modos_surf = sqrt(l(idx_modos_surf));
% q_modos_surf = sqrt((((b/a)*tan(k_modos_surf*L)).^2 + 1).*(k_modos_surf.^2));
% 
% % plot(q_modos_surf*2*Lx/pi, k_modos_surf, 'o', 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor','#D95319', 'DisplayName', 'Surface modes')
% plot([min(q) max(q)], [pi/(2*L) pi/(2*L)], 'Color', '#7E2F8E', 'HandleVisibility','off')%, 'DisplayName', '$\frac{\pi}{2 L}$')
% text(1,1.53,'$\frac{\pi}{2 L}$','FontSize', fsize, 'Color','#7E2F8E')
% % text(14,1.33,'Dispersion relation','FontSize', fsize, 'Color','#4169E1')
% % text(13.9,1.19,'Surface modes','FontSize', fsize, 'Color','#D95319')
% xlabel('$\frac{q L_x}{\pi}$', 'FontSize',fsize)
% ylabel('$k$', 'FontSize',fsize)
% %legend(ax, 'Location', 'northeastoutside')
% yticks([0 pi/8 pi/4 3*pi/8 pi/2])
% set(ax, 'yticklabels', {'0', '$\mathrm{\pi} /8$', '$\mathrm{\pi}/4$', '$3\mathrm{\pi}/8$', '$\mathrm{\pi}/2$'})
% ax.FontSize = fsize;
% xlim([0 Nmodos])
% grid


fsize = 28;
Nmodos = 20;
qn_interp = 0:1:Nmodos;
kn_interp = interp1(qn, k, qn_interp);

royalblue = '#0504aa';
firebrick = '#8f1402';
gray = '#3a3c3a';


figure;
ax = gca;
plot(qn_interp, kn_interp, '-o', 'color',royalblue, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', royalblue, 'DisplayName', 'Dispersion relation')
% ylim([min(k) pi/(2*L)])
xlim([0 Nmodos])
hold on

% Grafico los modos de superficie que salen de pdetool

idx_modos_surf = [9, 13, 17, 19, 21, 23, 26, 28, 31, 33, 34, 35, 36, 37, 39, 40, 41];
k_modos_surf = sqrt(l(idx_modos_surf));
q_modos_surf = sqrt((((b/a)*tan(k_modos_surf*L)).^2 + 1).*(k_modos_surf.^2));

plot(q_modos_surf*2*Lx/pi, k_modos_surf, 'o', 'MarkerFaceColor', firebrick, 'MarkerEdgeColor',firebrick, 'DisplayName', 'Surface modes')
plot([min(q) max(q)], [pi/(2*L) pi/(2*L)], '--','Color', gray, 'HandleVisibility','off')%, 'DisplayName', '$\frac{\pi}{2 L}$')
text(1,1.53,'$\frac{\pi}{2 L}$','FontSize', fsize, 'Color',gray)
text(15,1.36,'Dispersion relation','FontSize', fsize, 'Color',royalblue)
%text(14,1.33,'Dispersion relation','FontSize', fsize, 'Color',royalblue)
text(13.9,1.19,'Surface modes','FontSize', fsize, 'Color',firebrick)
xlabel('$\frac{q L_x}{\pi}$', 'FontSize',fsize+2)
ylabel('$k$', 'FontSize',fsize)
yticks([0 pi/8 pi/4 3*pi/8 pi/2])
set(ax, 'yticklabels', {'0', '$\mathrm{\pi} /8$', '$\mathrm{\pi}/4$', '$3\mathrm{\pi}/8$', '$\mathrm{\pi}/2$'})
ax.FontSize = fsize;
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
grid

fsize = 28;
Nmodos = 20;
qn_interp = 0:1:Nmodos;
kn_interp = interp1(qn, k, qn_interp);

royalblue = '#0504aa';
firebrick = '#8f1402';
gray = '#3a3c3a';


figure;
ax = gca;
plot(qn_interp, kn_interp, '-o', 'color',royalblue, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', royalblue, 'DisplayName', 'Dispersion relation')
% ylim([min(k) pi/(2*L)])
xlim([0 Nmodos])
hold on

% Grafico los modos de superficie que salen de pdetool

idx_modos_surf = [9, 13, 17, 19, 21, 23, 26, 28, 31, 33, 34, 35, 36, 37, 39, 40, 41];
k_modos_surf = sqrt(l(idx_modos_surf));
q_modos_surf = sqrt((((b/a)*tan(k_modos_surf*L)).^2 + 1).*(k_modos_surf.^2));

%plot(q_modos_surf*2*Lx/pi, k_modos_surf, 'o', 'MarkerFaceColor', firebrick, 'MarkerEdgeColor',firebrick, 'DisplayName', 'Surface modes')
plot([min(q) max(q)], [pi/(2*L) pi/(2*L)], '--','Color', gray, 'HandleVisibility','off')%, 'DisplayName', '$\frac{\pi}{2 L}$')
text(1,1.53,'$\frac{\pi}{2 L}$','FontSize', fsize, 'Color',gray)
text(15,1.36,'Dispersion relation','FontSize', fsize, 'Color',royalblue)
%text(14,1.33,'Dispersion relation','FontSize', fsize, 'Color',royalblue)
%text(13.9,1.19,'Surface modes','FontSize', fsize, 'Color',firebrick)
xlabel('$\frac{q L_x}{\pi}$', 'FontSize',fsize+2)
ylabel('$k$', 'FontSize',fsize)
yticks([0 pi/8 pi/4 3*pi/8 pi/2])
set(ax, 'yticklabels', {'0', '$\mathrm{\pi} /8$', '$\mathrm{\pi}/4$', '$3\mathrm{\pi}/8$', '$\mathrm{\pi}/2$'})
ax.FontSize = fsize;
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
grid


