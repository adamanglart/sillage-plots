import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from scipy.optimize import fsolve, root
from tqdm import tqdm
from scipy.signal import medfilt
from scipy.io import savemat
from figure_styling import *
from matplotlib import cm
plt.close('all')
plt.ion()

## Calculos siguiendo las cuentas de Piotr

g = 9.82 # m/s^2

def ky_deep_water(U,kx):
    # ecuacion 2.3
    return np.sqrt((U**4/g**2)*kx**4 - kx**2)


U = 0.4
kxmin = g/U**2
#kx = np.linspace(kxmin,200,1000)
kx = np.linspace(0,200,1000)
ky = ky_deep_water(U, kx)


def group_velocity(kx, ky):
    vx = np.sqrt(g)*kx/(2*(kx**2+ky**2)**(3/4))
    vy = np.sqrt(g)*ky/(2*(kx**2+ky**2)**(3/4))
    return vx, vy

vx, vy = group_velocity(kx, ky)

kxnorm = kx*U**2/g
kynorm = ky*U**2/g


## Angulo maximo
alpha_max = np.arcsin(1/3)
a = (1/4)*np.sin(alpha_max)
b = (1/4)*np.cos(alpha_max)

color1 = 'k'
color2 = 'firebrick'
color3 = 'lightcoral'

fig, ax = plt.subplots(figsize=(5,5))
ax.plot(kxnorm, kynorm  , c = color1)
ax.plot(kxnorm, -kynorm , c = color1)
ax.plot(-kxnorm, kynorm , c = color1)
ax.plot(-kxnorm, -kynorm, c = color1)
# j = 100
for j in range(100, len(kx),100):
    ax.arrow(kxnorm[j], kynorm[j] , vx[j]/U-1,  vy[j]/U,length_includes_head=True, head_width = 0.1, fill=True, color = color2)
    ax.arrow(kxnorm[j], -kynorm[j] , vx[j]/U-1,  -vy[j]/U,length_includes_head=True, head_width = 0.1, fill=True, color = color2)
ax.set_xlim((-6,6))
ax.set_ylim((-6,6))
ax.set_aspect('equal')
ax.set_xlabel(r'$\frac{k_x U^2}{g}$', fontsize = 14)
ax.set_ylabel(r'$\frac{k_y U^2}{g}$', fontsize = 14)
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/reldisp_sillage_deep.eps', dpi=100)


fig, ax = plt.subplots(figsize=(5,5))
polygon1 = Polygon([(1,0), (0.25,0), (0.25 + a ,b),], color='lightgray')
ax.add_patch(polygon1)
circle = plt.Circle(( 0.25 , 0 ), 0.25, fill=False, color='gray')
ax.add_artist(circle)
for i in range(0, len(vx), 150):
    ax.arrow(0,0, vx[i]/U,  vy[i]/U,length_includes_head=True, head_width = 0.01, fill=True, color = color3)
    ax.arrow(0,0, vx[i]/U, -vy[i]/U,length_includes_head=True, head_width = 0.01, fill=True, color = color3)
    ax.arrow(1,0, vx[i]/U-1,  vy[i]/U,length_includes_head=True, head_width = 0.01, fill=True, color = color2)
    ax.arrow(1,0, vx[i]/U-1, -vy[i]/U,length_includes_head=True, head_width = 0.01, fill=True, color = color2)

ax.text(0.1,0.25, '$V_g$', color=color3, fontsize=14)
ax.text(0.6,0.16, '$V_g \'$', color=color2, fontsize=14)
ax.set_xlim((0,1))
ax.set_ylim((-0.5,0.5))
ax.set_xticks([0,0.25,0.5,0.75,1])
ax.set_yticks([-0.5,-0.25,0,0.25,0.5])
ax.set_aspect('equal')
ax.set_xlabel(r'$V_x/U$', fontsize = 14)
ax.set_ylabel(r'$V_y/U$', fontsize = 14)
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/vgroup_deep.eps', dpi=100)



# PROFUNDIDAD FINITA



# def reldisp_profundidad_finita(k, Fr):
#     # ecuacion 2.9
#     norm_k = np.sqrt(k[0]**2 + k[1]**2)
#     return np.sqrt(norm_k*np.tanh(norm_k/(Fr**2))) - k[0]

def reldisp_profundidad_finita(ky, kx, Fr):
    # ecuacion 2.9
    norm_k = np.sqrt(kx**2 + ky**2)
    return np.sqrt(norm_k*np.tanh(norm_k/(Fr**2))) - kx


Frs = [0.3, 0.5, 0.6,  0.7, 1.0, 1.5, 2.0]

cmap  = cm.get_cmap('tab20b')
#cmap  = cm.get_cmap('Blues')
lines = ["-","--","-.",":", "-","--","-.",":"]
lines = [":","-.","--","-",":","-.","--","-"]

fig, ax = plt.subplots()
kynorm_fr = np.zeros((len(Frs), len(kxnorm)), dtype=float)
for j in range(len(Frs)):
    Fr = Frs[j]
    for i in range(len(kxnorm)):
        kynorm_fr[j,i] = fsolve(reldisp_profundidad_finita, x0 = kxnorm[i], args = (kxnorm[i], Fr))
    if Fr<1:
        kynorm_fr[j,kxnorm<0.99] = np.nan
    ax.plot(kxnorm,   kynorm_fr[j,:], c = cmap((j+1)/(len(Frs)+1)), label=str(round(Frs[j],2)))
    ax.plot(kxnorm,  -kynorm_fr[j,:], c = cmap((j+1)/(len(Frs)+1)))
    ax.plot(-kxnorm,  kynorm_fr[j,:], c = cmap((j+1)/(len(Frs)+1)))
    ax.plot(-kxnorm, -kynorm_fr[j,:], c = cmap((j+1)/(len(Frs)+1)))



# for j in range(len(Frs)):
#     if Fr<1:
#         kynorm_fr[j,kxnorm<0.99] = np.nan
#     ax.plot(kxnorm,   kynorm_fr[j,:], c = cmap(j/len(Frs)), label=str(Frs[j]))
#     ax.plot(kxnorm,  -kynorm_fr[j,:], c = cmap(j/len(Frs)))
#     ax.plot(-kxnorm,  kynorm_fr[j,:], c = cmap(j/len(Frs)))
#     ax.plot(-kxnorm, -kynorm_fr[j,:], c = cmap(j/len(Frs)))

color4 = 'firebrick'

ax.plot(kxnorm,   kynorm, '--', c = color4, label='deep water')
ax.plot(kxnorm,  -kynorm , '--', c = color4)
ax.plot(-kxnorm,  kynorm , '--', c = color4)
ax.plot(-kxnorm, -kynorm, '--', c = color4)
ax.set_xlim((-6,6))
ax.set_ylim((-6,6))
ax.set_aspect('equal')
ax.set_xlabel(r'$\frac{k_x U^2}{g}$', fontsize = 14)
ax.set_ylabel(r'$\frac{k_y U^2}{g}$', fontsize = 14)
plt.grid(True, alpha=0.6)
legend = ax.legend(bbox_to_anchor=(0.99, 0.95), edgecolor='None', facecolor='None')
ax.text(6.7,5.4, 'Fr', fontsize=14)
#legend.set_title("Fr")   
#plt.legend()
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/reldisp_sillage_hfinita.eps', dpi=100)



# Velocidad de grupo

def vx(kx, ky, Fr):
    norm_k = np.sqrt(kx**2 + ky**2)
    numerador   = (kx/norm_k)*np.tanh(norm_k/(Fr**2)) + (1/(Fr**2))*kx*(1/np.cosh(norm_k/(Fr**2)))**2
    denominador = 2*np.sqrt(norm_k*np.tanh(norm_k/(Fr**2)))
    return numerador/denominador - 1

def vy(kx, ky, Fr):
    norm_k = np.sqrt(kx**2 + ky**2)
    numerador   = (ky/norm_k)*np.tanh(norm_k/(Fr**2)) + (1/(Fr**2))*ky*(1/np.cosh(norm_k/(Fr**2)))**2
    denominador = 2*np.sqrt(norm_k*np.tanh(norm_k/(Fr**2)))
    return numerador/denominador



#Frs_for = np.linspace(0.001,4,1000)
Frs_for = np.arange(0.01,4,0.01)
#Frs_for = np.linspace(0.8,1.2,100)
#Frs_for = np.arange(0,1,0.1)
#Frs_for = [0.6,0.85,1.0,1.2,1.4]


kx_lim = 10
kx_inicial = np.linspace(0,kx_lim,3000)
kx_final = np.linspace(kx_lim,200,1000)
kx = np.concatenate((kx_inicial, kx_final))

ky_fr = np.zeros((len(Frs_for), len(kx)), dtype=float)
for j in tqdm(range(len(Frs_for))):
    Fr = Frs_for[j]
    for i in range(len(kx)):
        ky_fr[j,i] = fsolve(reldisp_profundidad_finita, x0 = kx[i], args = (kx[i], Fr))


cmap  = cm.get_cmap('tab20b')
alpha_max = np.zeros(len(Frs_for), dtype=float)
theta     = np.zeros(len(Frs_for), dtype=float)
#fig, ax = plt.subplots()
for j in range(len(Frs_for)):
    Fr = Frs_for[j]
    vx_fr = vx(kx, ky_fr[j,:], Fr)+1
    vy_fr = vy(kx, ky_fr[j,:], Fr)
    #ax.plot(vx_fr, vy_fr, '-',label=str(Fr))
    rs = vy_fr/(1-vx_fr)
    #idx_max = np.argmax(rs[~np.isnan(rs)])
    idx_max = np.argmax(medfilt(rs[~np.isnan(rs)],5))
    if Fr==1:
        print(idx_max)
    if Fr>=1:
        idx_max = idx_max +1
    #ax.plot(vx_fr[idx_max], vy_fr[idx_max], 'k.')
    alpha_max[j] = np.arctan(rs[idx_max])*180/np.pi
    theta[j]     = np.arctan(vy_fr[idx_max]/vx_fr[idx_max])*180/np.pi

# ax.set_aspect('equal')
# ax.set_xlim(xmax=1)
# plt.grid()
# plt.legend()
# plt.show()

Frs_total   = Frs_for.copy()
alpha_total = alpha_max.copy()
theta_total = theta.copy()

# plt.figure()
# plt.plot(Frs_total, alpha_total, 'k-')
# plt.xlabel('Fr')
# plt.ylabel(r'$\alpha$')
# plt.grid()
# plt.show()
# 
# plt.figure()
# plt.plot(Frs_total, theta, 'k-')
# plt.xlabel('Fr')
# plt.ylabel(r'$\theta$')
# plt.grid()
# plt.show()


plt.figure(figsize=(7,5))
plt.plot(Frs_total, alpha_total, '-', c='royalblue', label=r'$\alpha$')
plt.plot(Frs_total, theta_total, '-', c='firebrick', label=r'$\theta$')
plt.text(0.05, 15, 'Kelvin wake', color='royalblue')
plt.xlabel('Fr')
plt.ylabel('angle [degrees]')
plt.xlim(Frs_total[0], Frs_total[-1])
plt.legend(edgecolor='None', facecolor='None')
plt.grid(True, alpha=0.6)
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/angulos_sillage.eps', dpi=100)
np.save('/media/samantha/My Passport/these_sillage_plots/angulos_wake/alpha', alpha_total)
np.save('/media/samantha/My Passport/these_sillage_plots/angulos_wake/theta', theta_total)



stop
## Some Frs to illustrate how to obtain alpha

Frs_for = np.arange(0.6,1.6,0.2)

kx_lim = 10
kx_inicial = np.linspace(0,kx_lim,3000)
kx_final = np.linspace(kx_lim,200,1000)
kx = np.concatenate((kx_inicial, kx_final))

ky_fr = np.zeros((len(Frs_for), len(kx)), dtype=float)
for j in tqdm(range(len(Frs_for))):
    Fr = Frs_for[j]
    for i in range(len(kx)):
        ky_fr[j,i] = fsolve(reldisp_profundidad_finita, x0 = kx[i], args = (kx[i], Fr))

alpha_max = np.zeros(len(Frs_for), dtype=float)
fig, ax = plt.subplots()
for j in range(len(Frs_for)):
    Fr = Frs_for[j]
    vx_fr = vx(kx, ky_fr[j,:], Fr)+1
    vy_fr = vy(kx, ky_fr[j,:], Fr)
    if Fr<1:
        last_idx = np.argwhere(vy_fr>0.01)[0][0] -1
    else:
        last_idx = 0
    ax.plot(vx_fr[last_idx:], vy_fr[last_idx:], '-',label=str(round(Fr, 1)), c = cmap((j+1)/(len(Frs_for)+1)))
    rs = vy_fr/(1-vx_fr)
    idx_max = np.argmax(rs[~np.isnan(rs)])
    idx_max = np.argmax(medfilt(rs[~np.isnan(rs)],5))
    if Fr>=1:
        idx_max = idx_max +1
    #ax.plot(vx_fr[idx_max], vy_fr[idx_max], 'k.')
    ax.arrow(1,0, vx_fr[idx_max]-1, vy_fr[idx_max],length_includes_head=True, head_width = 0.01, fill=True, color = cmap((j+1)/(len(Frs_for)+1)))
    ax.arrow(0,0, vx_fr[idx_max], vy_fr[idx_max],length_includes_head=True, head_width = 0.01, fill=True, color = cmap((j+1)/(len(Frs_for)+1)), linestyle = '-', alpha = 0.6)
    alpha_max[j] = np.arctan(rs[idx_max])*180/np.pi
ax.set_aspect('equal')
ax.set_xlim(xmin=0, xmax=1)
ax.set_xlabel(r'$V_x/U$', fontsize = 14)
ax.set_ylabel(r'$V_y/U$', fontsize = 14)
legend = ax.legend(bbox_to_anchor=(0.99, 0.85), edgecolor='None', facecolor='None')
ax.text(1.08, 0.45, 'Fr', fontsize=14)
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/angulos_froude.eps', dpi=100)

