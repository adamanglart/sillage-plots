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
plt.close('all')
plt.ion()

path = '/media/samantha/My Passport/these_sillage_plots/mediciones_wake/'
mediciones_folder1 = 'mediciones-2023-06-06/'
mediciones_folder2 = 'mediciones-2023-06-22/'

g = 9.82    # m/s^2
h = 0.02    # m

#vs = np.arange(870,1010,10)
#vs = [880];
#vs = [900];
#vs = [1000];
#vs = [1100];

folder = [mediciones_folder2, mediciones_folder2, mediciones_folder1, mediciones_folder1]
vs = [880, 900, 1000, 1100]

letraA = ['(a)', '(c)', '(e)', '(g)']
letraB = ['(b)', '(d)', '(f)', '(h)']


angulos = [26.66, 26.58, 23.55, 20.88]


ii=3

#for ii in range(len(vs)):

mediciones_folder = folder[ii]
letra1 = letraA[ii]
letra2 = letraB[ii]

theta = angulos[ii]*np.pi/180



height_gra = np.load(path + mediciones_folder + ' gra_'+str(vs[ii])+'.npy')
height_mur = np.load(path + mediciones_folder + ' mur_'+str(vs[ii])+'.npy')

# print(height_mur.shape)

height_gra = height_gra[-1300:,:2034]
height_mur = height_mur[-1300:,:2034]


velocidad_dim = 0.36*vs[ii] - 2.55 # mm/s
froude_h      = (velocidad_dim/1000)/np.sqrt(g*h)

print(r'Fr = %.2f' % froude_h)
print('v = %.2f' % (velocidad_dim/10) +' cm/s')
print('v = %d' % vs[ii] +' ms/s')

ftp_proc_parameters = input_output.read_parameter_file(path+'processing_parameters.yaml')
pixel_size  = ftp_proc_parameters['MEASUREMENT']['pixel_size']
lin_min_idx = ftp_proc_parameters['MEASUREMENT']['lin_min_idx']
lin_max_idx = ftp_proc_parameters['MEASUREMENT']['lin_max_idx']
col_min_idx = ftp_proc_parameters['MEASUREMENT']['col_min_idx']
col_max_idx = ftp_proc_parameters['MEASUREMENT']['col_max_idx']


x = np.linspace(0, height_gra.shape[1]*pixel_size, height_gra.shape[1])
#y = np.linspace(0, height_gra.shape[0]*pixel_size, height_gra.shape[0])
y = np.linspace(0, -height_gra.shape[0]*pixel_size, height_gra.shape[0])


hmur = np.flipud(height_mur)
hgra = np.flipud(height_gra)


mask = np.zeros_like(hmur)
row = 617
col = 1837
radius = 62

rr, cc = disk((row, col), radius)
mask[rr, cc] = 1
mask[592:642, 1892:] = 1

grating = np.zeros((104, hmur.shape[1]))
mur     = np.ones((104, hmur.shape[1]))
idx_inicial_grating = 30
length_grating = 56
periodicity_grating = 2*length_grating
for i in range(hmur.shape[1]//periodicity_grating):
    grating[:-15, idx_inicial_grating + i*periodicity_grating : idx_inicial_grating + i*periodicity_grating + length_grating] = 1

hmur_extended = np.concatenate((mur, hmur))
hgra_extended = np.concatenate((grating, hgra))




mask_mur_extended = np.concatenate((mur, mask))
mask_gra_extended = np.concatenate((grating, mask))

x = np.linspace(0, hgra_extended.shape[1]*pixel_size, hgra_extended.shape[1])
#y = np.linspace(0, height_gra.shape[0]*pixel_size, height_gra.shape[0])
y = np.linspace(104*pixel_size, -(hgra_extended.shape[0]-104)*pixel_size, hgra_extended.shape[0])



fig, ax = plt.subplots(1,2, sharey=True, figsize=(10,4))#, dpi=200)

im0 = ax[0].imshow(ma.masked_array(hmur_extended, mask = mask_mur_extended), aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], y[0]])
im1 = ax[1].imshow(ma.masked_array(hgra_extended, mask = mask_mur_extended),  aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], y[0]])
# im0 = ax[0].imshow(hmur, aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], 0])
# im1 = ax[1].imshow(hgra,  aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], 0])

ax[0].axhline(y = 0, color='k')

A = 0.08
ax[0].arrow(0.32, -0.08, A*np.cos(theta), A*np.sin(theta), width = 0.0004, head_width=0.006, facecolor='k', edgecolor='k')
ax[0].arrow(0.12, -0.08 + A*np.sin(theta), A*np.cos(theta), -A*np.sin(theta), width = 0.0004, head_width=0.006, facecolor='k', edgecolor='k')
ax[0].text(0.14, -0.085, '$R$')
ax[0].text(0.36, -0.085, '$\mathbf{k}$')

ax[1].arrow(0.32, -0.08, A*np.cos(theta), A*np.sin(theta), width = 0.0004, head_width=0.006, facecolor='k', edgecolor='k')
ax[1].arrow(0.12, -0.08 + A*np.sin(theta), A*np.cos(theta), -A*np.sin(theta), width = 0.0004, head_width=0.006, facecolor='k', edgecolor='k')
ax[1].text(0.14, -0.085, '$R$')
ax[1].text(0.36, -0.085, '$\mathbf{k}$')



for i in range(hmur.shape[1]//periodicity_grating):
    rect = patches.Rectangle((pixel_size*(idx_inicial_grating + i*periodicity_grating), pixel_size*15), pixel_size*length_grating, pixel_size*89, linewidth=1, edgecolor='k', facecolor='k')
    ax[1].add_patch(rect)

ax[0].set_title('solid wall')
ax[1].set_title('grating')

ax[0].text(0.015, -0.35, letra1)
ax[1].text(0.015, -0.35, letra2)

im0.set_clim((-4.5,4.5))
im1.set_clim((-4.5,4.5))
ax[0].set_ylabel('y [m]')
ax[0].set_xlabel('x [m]')
ax[1].set_xlabel('x [m]')

cbar = fig.colorbar(im0, ax=ax.ravel().tolist(), fraction=0.046, pad=0.04, shrink = 0.72)
cbar.ax.set_ylabel('h [mm]', rotation=270, labelpad=20)

plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/sillage_experimental_' + '%d' % vs[ii]  +'.png')
plt.close()





#plt.savefig(path_graficos + '%d' % vs[ii], dpi=200)

# stop
# 
# fig, ax = plt.subplots(1,2, sharey=True, figsize=(8.27,4))#, dpi=200)
# 
# im0 = ax[0].imshow(ma.masked_array(hmur, mask = mask), aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], 0])
# im1 = ax[1].imshow(ma.masked_array(hgra, mask = mask),  aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], 0])
# # im0 = ax[0].imshow(hmur, aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], 0])
# # im1 = ax[1].imshow(hgra,  aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], 0])
# 
# 
# ax[0].set_title('solid wall')
# ax[1].set_title('grating')
#     
# im0.set_clim((-5,5))
# im1.set_clim((-5,5))
# ax[0].set_ylabel('y [m]')
# ax[0].set_xlabel('x [m]')
# ax[1].set_xlabel('x [m]')
# 
# cbar = fig.colorbar(im0, ax=ax.ravel().tolist(), fraction=0.046, pad=0.04, shrink = 0.56)
# # cbar.ax.set_ylabel('h [mm]', rotation=270)
# plt.show()
# #plt.savefig(path_graficos + '%d' % vs[ii], dpi=200)
# 
# 
# mask = np.zeros_like(hmur)
# row = 617
# col = 1837
# radius = 62
# 
# 
# rr, cc = disk((row, col), radius)
# mask[rr, cc] = 1
# mask[592:642, 1892:] = 1
# 
# plt.figure()
# plt.imshow(1-grating, cmap='gray')
# plt.show()
# 
# 
# 
# 
# 
# plt.figure()
# plt.imshow(ma.masked_array(hmur, mask = mask))
# #plt.imshow(mask)
# plt.show()
# 
# 
# 
# gray = sio.imread('/media/samantha/SP PHD U3/mediciones-2023-06-22/grating/gray/im_C001H001S0001000001.tif')
# plt.figure()
# plt.imshow(gray, aspect='auto', cmap='gray')
# plt.axhline(y = 1300)
# plt.clim(0,500)
# plt.grid()
# plt.show()
# 
