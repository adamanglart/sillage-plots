clear all
close all
clc


set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultLineLineWidth',1.5)

fsize = 16;


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

ks = linspace(0.001,pi,100)./a;
kx = ks.*cos(theta);
ky = ks.*sin(theta);


%ks = [pi/a];
%ks = [0.32];
%ks = [pi/2];

%fig= figure('units','inch','position',[10,10,10,7]);
bs = [0.2,0.5,0.7,1];
bs = 0.02:0.02:0.6;
%bs = [0.2];
for jj=1:length(bs)
b = bs(jj)


%Ls = 1:1:5;
% Ls = [1];
% subplot(2,2,jj);
% for ll=1:length(Ls)
% 
% L = Ls(ll);
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



reflection(bb) = R(N);

end

%% Fit of the complex resonance frequency

M = length(ks);
%[val, idx] = max(diff(unwrap(angle(reflection(1:floor(0.9*M)))))); % la resonancia ocurre donde la pendiente de la fase es maxima
idx = M/2;
kA0 = ks(idx);



clear difference
thetas = linspace(-pi,pi, 100);
alpha0 = b/cos(theta);
alphas = linspace(-0.5*alpha0, 0.5*alpha0, 200);

for ii=1:length(thetas)
	for kk=1:length(alphas)
		theta_i = thetas(ii);
		alphaR = alphas(kk);
		Rexpression = exp(1i*theta_i).*(ks - kA0 - 1i*alphaR)./(ks - kA0 + 1i*alphaR);
		difference(ii,kk) = norm(R-Rexpression);

	end
end



% if jj==10
% 	figure;
% 	contourf(alphas, thetas, difference, 20, 'LineStyle','none')
% 	ylabel('$\theta$')
% 	xlabel('$\alpha_R$')
% 	set(gca,'YTick',-pi:pi/2:pi) 
% 	set(gca,'YTickLabel',{'','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})
% 	set(gca,'YTickLabel',{'$-\pi$','$-\pi/2$','$0$', '$\pi/2$','$\pi$'})
% 	cmap = cmocean('ice');
% 	set(gca, 'colormap', cmap);
% 	shading interp
% 	colorbar('TickLabelInterpreter', 'latex');
% 	title(join(['b/a = ',num2str(b)]))
% 	set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
% end

[m n] = min(difference(:));
[idx_x idx_y] = ind2sub(size(difference), n);
aA0 = alphas(idx_y);
tA0 = thetas(idx_x);

kR(jj) = kA0;
aR(jj) = aA0;
tR(jj) = tA0;


W = 20;
idx_min = idx - W;
idx_max = idx + W;

if idx_min<1, idx_min=1; end
if idx_max>M, idx_max=M; end

kcut = ks(idx_min:idx_max);
Rcut = reflection(idx_min:idx_max);


%% specular reflection
num = 1j.*ky - ks.*b.*tan(ks.*L);
den = 1j.*ky + ks.*b.*tan(ks.*L);
Rspecular = num./den;

Rspecular_cut = Rspecular(idx_min:idx_max);



%% Fitteo con esta funcion (ec. 4 del articulo)
funRk=@(fit,k)((k-fit(1))./(k-conj(fit(1)))).*exp(1j*real(fit(2)));
% Modal
eta_kfitRk = lsqcurvefit(funRk, [pi/2 + 1j*b/cos(theta), tR(jj)], kcut, Rcut);
eta_kRk=funRk(eta_kfitRk,kcut);

kR_fit(jj) = real(eta_kfitRk(1));
aR_fit(jj) = imag(eta_kfitRk(1));
tR_fit(jj) = eta_kfitRk(2);

% Rspecular
% [val_s,idx_s]=min(abs(ks-kR_fit(jj)));
% W_s = 15;
% idx_min_s = idx_s - W_s;
% idx_max_s = idx_s + W_s;


% if jj<11
% 	idx_min_s = idx_min;
% 	idx_max_s = idx_max;
% else
% 	idx_min_s = idx_min_s - 1;
% 	idx_max_s = idx_min_s + 1;
% end

idx_min_s = 1;
idx_max_s = M-1;




kcut_spec = ks(idx_min_s:idx_max_s);
Rspecular_cut = Rspecular(idx_min_s:idx_max_s);

%eta_specular_kfitRk = lsqcurvefit(funRk, [eta_kfitRk(1), 0], kcut, Rspecular_cut);
%eta_specular_kfitRk = lsqcurvefit(funRk, [pi/2 + 1j*b/cos(theta), 0], kcut, Rspecular_cut);
%eta_specular_kfitRk = lsqcurvefit(funRk, [eta_kfitRk(1), 0], kcut_spec, Rspecular_cut);
%eta_specular_kfitRk = lsqcurvefit(funRk, [pi/2 + 1j*b/cos(theta), 0], kcut_spec, Rspecular_cut);
eta_specular_kfitRk = lsqcurvefit(funRk, [real(eta_kfitRk(1)) + 1j*b/cos(theta), 0], kcut_spec, Rspecular_cut);
%eta_specular_kfitRk = lsqcurvefit(funRk, [pi/2 + 1j*b/cos(theta), eta_kfitRk(2)], kcut_spec, Rspecular_cut);
eta_specular_kRk=funRk(eta_specular_kfitRk,kcut_spec);

kR_specular_fit(jj) = real(eta_specular_kfitRk(1));
aR_specular_fit(jj) = imag(eta_specular_kfitRk(1));
tR_specular_fit(jj) = eta_specular_kfitRk(2);




%if jj>15
if jj==10

	kcut_spec = ks(idx_min:idx_max);
	Rspecular_cut = Rspecular(idx_min:idx_max);
	eta_specular_kfitRk = lsqcurvefit(funRk, [real(eta_kfitRk(1)) + 1j*b/cos(theta), 0], kcut_spec, Rspecular_cut);
	eta_specular_kRk=funRk(eta_specular_kfitRk,kcut_spec);

	fig = figure('units','inch','position',[10,10,8,4]);
	plot(ks, unwrap(angle(reflection)), 'DisplayName', 'multimodal method', 'color', royalblue);%, 'color', CM(ll,:))
	hold on
	plot(kcut, unwrap(angle(eta_kRk)), 'k--', 'HandleVisibility','off')
	plot(ks, unwrap(angle(Rspecular)), 'DisplayName', 'effective boundary condition', 'color', firebrick);%, 'color', CM(ll,:))
	plot(kcut_spec, unwrap(angle(eta_specular_kRk)), 'k--','DisplayName', 'fit')
	set(gca,'XTick',0:pi/4:pi) 
	set(gca,'YTick',0:pi/4:2*pi) 
	xlabel('$ka$')
	set(gca,'XTickLabel',{'0', '$\pi/4$','$\pi/2$', '$3\pi/2$', '$\pi$'})
	ylabel('phase($R$)')
	set(gca,'YTickLabel',{'0', '$\pi/4$','$\pi/2$','$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$'})
	lgd = legend('Location','northwest', 'Box', 'off', 'edgecolor', 'white');
	%lgd.Title.String = '$b/a$';
	xlim([0 pi])
	ylim([0 2*pi])
	grid
	set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
	%title(join(['$b/a$ = ', num2str(b)]), 'interpreter','latex')
	saveas(fig,'/home/samantha/Dropbox/PhD/these/figures/wake/phaseR_fit_b02.png')
end

end

%end


bs_dense = 0:0.01:0.61;


fig = figure('units','inch','position',[10,10,10,7]);
plot(bs, kR_fit, 'o', 'color', royalblue, 'MarkerFaceColor', royalblue, 'DisplayName', 'multimodal method')
hold on
plot(bs, aR_fit, 'ro', 'color', firebrick, 'MarkerFaceColor', firebrick,  'HandleVisibility','off')
plot(bs, tR_fit, 'go', 'color', green, 'MarkerFaceColor', green, 'HandleVisibility','off' )

plot(bs, kR_specular_fit, 's', 'color', royalblue,  'MarkerFaceColor', 'white', 'DisplayName', 'effective boundary condition')
plot(bs, aR_specular_fit, 'rs', 'color', firebrick, 'MarkerFaceColor', 'white', 'HandleVisibility','off' )
plot(bs, tR_specular_fit, 'gs', 'color', green,     'MarkerFaceColor', 'white',  'HandleVisibility','off')

plot(bs_dense, (pi/2).*ones([1,length(bs_dense)]), '--', 'color', royalblue, 'DisplayName', 'zero-pole approximation')
plot(bs_dense, bs_dense/cos(theta), '--', 'color', firebrick,  'HandleVisibility','off')
plot(bs_dense, 0.*ones([1,length(bs_dense)]), '--', 'color', green,  'HandleVisibility','off')

xlabel('$b/a$')
ylabel('')
xlim([0 0.61])
lgd = legend('Location','northwest', 'Box', 'off', 'edgecolor', 'white');
lgd.Position(1) = 0.2;
lgd.Position(2) = 0.6;
inicio_texto = 0.62;
delta_texto = 0.12;
hTxt(1,1)=text(-0.05,inicio_texto,{'$k_R$,'},'Rotation',90,'color',royalblue);
hTxt(2,1)=text(-0.05,inicio_texto + delta_texto,{'$\alpha_r,$'},'Rotation',90,'color', firebrick);
hTxt(3,1)=text(-0.05,inicio_texto + 2*delta_texto,{'$\varphi$'},'Rotation',90,'color',green);
set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
hold off
grid
saveas(fig,'/home/samantha/Dropbox/PhD/these/figures/wake/kres_modal.png')

% fig = figure('units','inch','position',[10,10,10,7]);
% plot(bs, kR_fit, 'o', 'color', royalblue, 'MarkerFaceColor', royalblue)
% hold on
% plot(bs, aR_fit, 'ro', 'color', firebrick, 'MarkerFaceColor', firebrick)
% plot(bs, tR_fit, 'go', 'color', green, 'MarkerFaceColor', green)
% xlabel('$b/a$')
% ylabel('')
% inicio_texto = 0.62;
% delta_texto = 0.12;
% hTxt(1,1)=text(-0.05,inicio_texto,{'$k_R$,'},'Rotation',90,'color',royalblue);
% hTxt(2,1)=text(-0.05,inicio_texto + delta_texto,{'$\alpha_r,$'},'Rotation',90,'color', firebrick);
% hTxt(3,1)=text(-0.05,inicio_texto + 2*delta_texto,{'$\varphi$'},'Rotation',90,'color',green);
% set(findall(gcf,'-property','FontSize'),'FontSize',fsize)
% hold off
% grid
