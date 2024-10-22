import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from figure_styling import *
import sys
sys.path.insert(0, '/home/samantha/Documents/')
from FTP_PhD import input_output
from skimage.draw import disk
import skimage.io as sio
import matplotlib.patches as patches
from scipy.signal import find_peaks
from scipy.optimize import fsolve, root
plt.close('all')
plt.ion()

path = '/media/samantha/My Passport/these_sillage_plots/mediciones_wake/'
mediciones_folder = 'mediciones-2023-06-22/'

g = 9.82    # m/s^2
h = 0.02    # m

vs = np.arange(870,1010,10)
vs = [900];
ii=0

#for ii in range(len(vs)):


height_mur = np.load(path + mediciones_folder + ' mur_'+str(vs[ii])+'.npy')


velocidad_dim = 0.36*vs[ii] - 2.55 # mm/s
froude_h      = (velocidad_dim/1000)/np.sqrt(g*h)

print(r'Fr = %.2f' % froude_h)
print('v = %.2f' % (velocidad_dim/10) +' cm/s')
print('v = %d' % vs[ii] +' ms/s')

ftp_proc_parameters = input_output.read_parameter_file(path+'processing_parameters.yaml')
pixel_size  = ftp_proc_parameters['MEASUREMENT']['pixel_size']
# lin_min_idx = ftp_proc_parameters['MEASUREMENT']['lin_min_idx']
# lin_max_idx = ftp_proc_parameters['MEASUREMENT']['lin_max_idx']
# col_min_idx = ftp_proc_parameters['MEASUREMENT']['col_min_idx']
# col_max_idx = ftp_proc_parameters['MEASUREMENT']['col_max_idx']

hmur = np.flipud(height_mur)
height  = hmur[:,:1750]

x = np.linspace(0,  height.shape[1]*pixel_size, height.shape[1])
y = np.linspace(0, -height.shape[0]*pixel_size, height.shape[0])



# fig, ax = plt.subplots()
# im = ax.imshow(height,aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], y[0]])
# im.set_clim((-4.5,4.5))
# ax.set_ylabel('$y$ [m]')
# ax.set_xlabel('$x$ [m]')
# cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.63)
# #cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.79)
# cbar.ax.set_ylabel('h [mm]', rotation=270, labelpad=20)
# 
# plt.tight_layout()
# plt.show()


def ky_deep_water(U,kx):
    # ecuacion 2.3
    return np.sqrt((U**4/g**2)*kx**4 - kx**2)



def reldisp_profundidad_finita_dim(ky, kx, U, h):
    # ecuacion 4.35
    norm_k = np.sqrt(kx**2 + ky**2)
    return np.sqrt(g*norm_k*np.tanh(norm_k*h)) - U*kx

fourier = np.fft.fftshift(np.fft.fftshift(np.fft.fft2(height), axes=0), axes=1)
kx = np.fft.fftshift(np.fft.fftfreq(len(x), np.diff(x)[0]))
ky = np.fft.fftshift(np.fft.fftfreq(len(y), np.diff(y)[0]))
ft = np.log(np.abs(fourier))/np.max(np.log(np.abs(fourier)))


U = velocidad_dim/1000
Fr = froude_h


kx_reldisp = np.linspace(0, np.max(kx), 2000)
ky_deep = ky_deep_water(U, kx_reldisp)



ky_reldisp_finita = np.zeros(len(kx_reldisp), dtype=float)
for i in range(len(kx_reldisp)):
    ky_reldisp_finita[i] = fsolve(reldisp_profundidad_finita_dim, x0 = kx_reldisp[i], args = (kx_reldisp[i], U, h))


fig, ax = plt.subplots()
im = ax.imshow(ft, aspect='auto', cmap=cmocean.cm.ice_r, extent = [kx[0], kx[-1], ky[0], ky[-1]])

ax.plot( kx_reldisp,   ky_reldisp_finita, '--', c ='k')
ax.plot( kx_reldisp,  -ky_reldisp_finita, '--', c ='k')
ax.plot(-kx_reldisp,   ky_reldisp_finita, '--', c ='k')
ax.plot(-kx_reldisp,  -ky_reldisp_finita, '--', c ='k')


ax.plot( kx_reldisp,   ky_deep, '--', c ='k')
ax.plot( kx_reldisp,  -ky_deep, '--', c ='k')
ax.plot(-kx_reldisp,   ky_deep, '--', c ='k')
ax.plot(-kx_reldisp,  -ky_deep, '--', c ='k')


ax.set_xlim([-150, 150])
ax.set_ylim([-150, 150])
im.set_clim([0,1])
cbar = fig.colorbar(im, ax=ax)#, fraction=0.046, pad=0.04, shrink = 0.63)
plt.tight_layout()
plt.show()








