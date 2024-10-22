clear all
close all
clc


set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultLineLineWidth',1.5)

fsize = 14;



% domain size
l1 = 5;

% grating size
a = 1;
b = 0.2;
L = 1; 
l2 = L;


% incident wave
%k = 2.18;
theta = 35*pi/180;

ks = linspace(0.001,pi,100)./a;
kx = ks.*cos(theta);
ky = ks.*sin(theta);


ks = [(9/10)*pi/a];
%ks = [0.32];
%ks = [pi/2];

%fig= figure('units','inch','position',[10,10,10,7]);
%bs = [0.2,0.5,0.7,1];
bs = [0.5];
for jj=1:length(bs)
b = bs(jj)


%Ls = 1:1:5;
Ls = [5];
%subplot(2,2,jj);
for ll=1:length(Ls)

L = Ls(ll);
ls = L;
for bb=1:length(ks)

k = ks(bb);

alpha = k*sin(theta); % kx0
beta = k*cos(theta); %ky0

kinc2 = alpha^2 + beta^2;


% number of modes
%N1 = 29;
%N1 = 21;
N1 = 5;
N2 = N1;
N = N1;
%N2 = 2*N1 - 1;
%N2 = floor(N1/2);


% modes in x
%bn2 = (0:N1-1)*pi/b;
bn1 = alpha + 2*(0:N1-1)*pi/a;
% bn2 = (-(N2-1):N2-1)*pi/b;

%bn1 = alpha + 2*(-(N1-1):N1-1)*pi/a;
bn2 = (0:N2-1)*pi/b;

bn1 = alpha + 2*(-(N1-1):N1-1)*pi/a;
bn2 = (0:2*N2-1)*pi/b;


% modes in y
kn1 = sqrt(kinc2 -bn1.^2);
kn2 = sqrt(kinc2 -bn2.^2);


N1  = length(kn1);
N2  = length(kn2);

norm1(1:N1) = 1/sqrt(a); 
norm2(1) = 1/sqrt(b); 
norm2(2:N2) = sqrt(2/b);


x=linspace(0,b,1e3);
for nn=1:N1,
	for mm=1:N2, 
        	phip=norm1(nn)*exp(-1j*bn1(nn).*x); 
		phim=norm2(mm)*cos(bn2(mm).*x);
        	C(mm,nn)=trapz(x,phip.*phim);
    end
end

sizeC = size(C);


% dividiendo por cos(kL)
% H = [-conj(C), eye(sizeC(1));
% diag(1i*kn1), kn2.*tan(kn2*L).*C'];
H = [-conj(C), eye(sizeC(1));
diag(1i*kn1),  C.'*diag(kn2.*tan(kn2*L))];

% H = [-conj(C), cos(kn2*L).*eye(sizeC(1));
% diag(1i*kn1),  C.'*diag(kn2.*sin(kn2*L))];
% %diag(1i*kn1), kn2.*sin(kn2*L).*C'];


S1 = conj(C(:,N));
%S1 = conj(C(:,1));
S2 = zeros(size(kn1))';
%S2(1) = 1i*beta;
S2(N) = 1i*beta;
V = [S1; S2];

sol = H\V;

R = sol(1:N1);
A = sol(N1+1:end);

% Region 1
Nx1=400; Ny1=100;
x1=ones(Ny1,1)*linspace(0,a,Nx1);
y1=linspace(-l1,0,Ny1)'*ones(1,Nx1);

phi1 = norm1(1)*exp(1i*beta*y1).*exp(1j*alpha*x1);
for ii=1:N1,   
    phi1 = phi1 + R(ii)*exp(-1i*kn1(ii)*y1)*norm1(ii).*exp(1i*bn1(ii)*x1);   
end

% Region 2
Nx2=400; Ny2=100;
x2=ones(Ny2,1)*linspace(0,b,Nx2);
y2=linspace(0,l2,Ny2)'*ones(1,Nx2);

phi2 = zeros(size(y2));
for ii=1:N2,   
    % dividiendo por cos(kL)
    phi2 = phi2 + A(ii)*(cos(kn2(ii)*(y2-L))/cos(kn2(ii)*L))*norm2(ii).*cos(bn2(ii)*x2);  
    % figure;
    % pcolor(x2, y2, real(phi2))
    % title(num2str(ii))
    %phi2 = phi2 + A(ii)*(cos(kn2(ii)*(y2-L)))*norm2(ii).*cos(bn2(ii)*x2);   
end   

% figure;
% pcolor(x1,y1,real(phi1))
% hold on 
% pcolor(x2,y2,real(phi2)) 
% axis auto , colormap(cmocean('balance')), shading interp, axis equal, colorbar
% title('Re($\phi$)','fontsize',20)
% % 
%% Periodic structure

Nperiods = 10;
fig= figure('units','inch','position',[10,10,8.5,5]);
pcolor(x1,y1,real(phi1))
hold on 
pcolor(x2,y2,real(phi2))
for kk=1:Nperiods
	pcolor(x1+kk*a,y1,real(phi1*exp(1i*alpha*a*kk)))
	pcolor(x2+kk*a,y2,real(phi2*exp(1i*alpha*a*kk))) 
end
%title('Re($\phi$)','fontsize',20)
xlabel('$x/a$')
ylabel('$y/a$')
axis auto , colormap(cmocean('balance')), shading interp, axis equal, colorbar
%title('Re($\phi$)','fontsize',20)
ax = gca;
ax.FontSize = fsize;
set(ax, 'CLim', [-2 2])
%set(ax, 'CLim', [-1.5 1.5])
cbh = colorbar(ax, 'TickLabelInterpreter', 'latex');
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
ylim([-l1 l2])
% xlim([-l1 l2])
% xticks(-5:2.5:5)
% ylim([0 Nperiods]) 
grid('off')
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[7.5, 5]);

x1p = 2.82;
y1p = -3.53;
%x2p = 3.82;
%y2p = -1;
x2p = x1p + 3*sin(theta);
y2p = y1p + 3*cos(theta);
arh2 = annotation('arrow', 'HeadStyle', 'plain', 'LineWidth', 1,'HeadLength', 5, 'HeadWidth', 5, 'color', 'k');
arh2.Units = 'normalized';
arh2.Parent = fig.CurrentAxes;
arh2.Position = [x1p, y1p, x2p-x1p, y2p-y1p];
%arh2.Position = [0.3, 0.3, sin(theta), cos(theta)];
arh2.Color = 'k';

text(3.42, -2.26, '\textbf{k}', 'FontSize', 14)






%print(fig,'/home/samantha/Dropbox/PhD/these/figures/wake/field_modal','-dpdf','-r300')
saveas(fig,'/home/samantha/Dropbox/PhD/these/figures/wake/field_modal.png')



reflection(bb) = R(N);

end

% CM = autumn(length(Ls));
% %CM = cmocean('balance')(length(Ls));
% 
% %% specular reflection
% num = 1j.*ky - ks.*b.*tan(ks.*L);
% den = 1j.*ky + ks.*b.*tan(ks.*L);
% Rspecular = num./den;
% 
% 
% %figure;
% plot(ks, unwrap(angle(reflection)), 'DisplayName', join(['L/a = ', num2str(L)]), 'color', CM(ll,:))
% hold on
% plot(ks, unwrap(angle(Rspecular)), '--', 'HandleVisibility','off', 'color', CM(ll,:))
% %plot(ks, unwrap(angle(Rspecular)), '--', 'color', CM(ll,:))
% if jj==2
% 	legend('Location','southeast')
% end
% set(gca,'XTick',0:pi/4:pi) 
% set(gca,'YTick',0:pi/4:2*pi) 
% set(gca,'XTickLabel',{'', '','', '', ''})
% set(gca,'YTickLabel',{'', '','','', '', '', '', '', ''})
% if jj==3 | jj==4
% 	xlabel('$ka$')
% 	set(gca,'XTickLabel',{'0', '$\pi/4$','$\pi/2$', '$3\pi/2$', '$\pi$'})
% end
% if jj==1 | jj==3
% 	ylabel('phase($R$)')
% 	set(gca,'YTickLabel',{'0', '$\pi/4$','$\pi/2$','$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$'})
% end
% 
% xlim([0 pi])
% ylim([0 2*pi])
% grid
% set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
% %title('Modal method', 'interpreter','latex')
% title(join(['b/a = ', num2str(b)]), 'interpreter','latex')

end

end

% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[10, 7]);
% print(fig,'/home/samantha/Dropbox/PhD/these/figures/wake/comparacion_modal_bc_efectiva','-dpdf','-r300')


% figure;
% pcolor(x1,y1,imag(phi1))
% hold on 
% pcolor(x2,y2,imag(phi2)) 
% axis auto , colormap('jet'), shading interp, axis equal, colorbar
% title('Im($\phi$)','fontsize',20)
% 
% 
% figure;
% pcolor(x1,y1,abs(phi1))
% hold on 
% pcolor(x2,y2,abs(phi2)) 
% axis auto , colormap('jet'), shading interp, axis equal, colorbar
% title('$| \phi |$','fontsize',20)


