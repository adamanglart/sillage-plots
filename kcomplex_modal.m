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


royalblue = '#0504aa';
firebrick = '#8f1402';
gray = '#3a3c3a';
green = '#006400';

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

ks_real = linspace(0.001,pi,100)./a;
ks_imag = linspace(-0.5*pi, 0.5*pi, 100);
% ks = ks_real + 1j.*ks_imag;
% kx = ks.*cos(theta);
% ky = ks.*sin(theta);


%ks = [pi/a];
%ks = [0.32];
%ks = [pi/2];

%fig= figure('units','inch','position',[10,10,10,7]);
bs = [0.2,0.5,0.7,1];
%bs = 0.02:0.02:0.6;
bs = [0.3];
%bs = [0.2, 0.3, 0.4];
for jj=1:length(bs)
b = bs(jj)


%Ls = 1:1:5;
% Ls = [1];
% subplot(2,2,jj);
% for ll=1:length(Ls)
% 
% L = Ls(ll);
ls = L;
for bb=1:length(ks_real)

for dd=1:length(ks_imag)

k = ks_real(bb) + 1j*ks_imag(dd);

alpha = k*sin(theta); % kx0
beta = k*cos(theta); %ky0

kinc2 = alpha^2 + beta^2;


% number of modes
%N1 = 29;
%N1 = 21;
N1 = 5;
%N1 = 5;
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



%kn1 = kn1(imag(kn1)<0);
idx_not_ok = find(imag(kn1)>0);
kn1(idx_not_ok) = conj(kn1(idx_not_ok));





%stop


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
% % Periodic structure
% Nperiods = 10;
% figure;
% pcolor(x1,y1,real(phi1))
% hold on 
% pcolor(x2,y2,real(phi2))
% for kk=1:Nperiods
% 	pcolor(x1+kk*a,y1,real(phi1*exp(1i*alpha*a*kk)))
% 	pcolor(x2+kk*a,y2,real(phi2*exp(1i*alpha*a*kk))) 
% end
% axis auto , colormap(cmocean('balance')), shading interp, axis equal, colorbar
% %title('Re($\phi$)','fontsize',20)



reflection(bb, dd) = R(N);


end

bb
end

%[kR, kI] = meshgrid(ks_real,ks_imag);

figure;
pcolor(ks_real, ks_imag, 20*log(abs(reflection)).')
hold on
%5contourf(ks_real, ks_imag, 20*log(abs(reflection)).', 'k', 'levels', 4)
%plot(ks_real, 0*ks_real, 'k--')
hold off
%dx = 0.2
%xlim([1.46-dx, 1.46+dx])
%ylim([-1, 1])
cmap = cmocean('balance');
set(gca, 'colormap', cmap);
shading flat
set(gca, 'clim', [-50 50]);
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})
set(gca,'YTick',-pi/2:pi/4:pi/2) 
set(gca,'YTickLabel',{'$-\pi/2$','$-\pi/4$', '$0$', '$\pi/4$','$\pi/2$'})
xlabel('Re($ka$)')
ylabel('Im($ka$)')
title(join(['b/a = ', num2str(round(b, 2))]))
colorbar('TickLabelInterpreter', 'latex');
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
%savefig('/media/samantha/SP PHD U3/simulaciones-2023-06-22/R05.png')
%save(join(['/media/samantha/My Passport/these_sillage_plots/R_kcomplex/', num2str(round(b, 2)), '.mat']), 'reflection')

jj
end
