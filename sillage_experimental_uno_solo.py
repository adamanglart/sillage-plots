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
plt.close('all')
plt.ion()

path = '/media/samantha/My Passport/these_sillage_plots/mediciones_wake/'
mediciones_folder = 'mediciones-2023-06-22/'

g = 9.82    # m/s^2
h = 0.02    # m

vs = np.arange(870,1010,10)
vs = [900];
#vs = [880];
ii=0

#for ii in range(len(vs)):


height_gra = np.load(path + mediciones_folder + ' gra_'+str(vs[ii])+'.npy')
height_mur = np.load(path + mediciones_folder + ' mur_'+str(vs[ii])+'.npy')


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



fig, ax = plt.subplots()
im = ax.imshow(ma.masked_array(hmur, mask=mask),aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], y[0]])
im.set_clim((-4.5,4.5))
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.63)
#cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.79)
cbar.ax.set_ylabel('h [mm]', rotation=270, labelpad=20)

plt.tight_layout()
plt.show()
# plt.savefig('/home/samantha/Desktop/fig1.png', dpi=200)



# fig, ax = plt.subplots()
# im = ax.imshow(ma.masked_array(hmur, mask=mask),aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], y[0]])
# #ax.axhline(y = -0.175)
y0 = y[615]
# ax.axhline(y = y0, color='k')
m = np.tan(26.81*np.pi/180)
# idx_x0 = 350 # 0.1
idx_x0 = 244
x0 = x[idx_x0]
# ax.plot(x, m*(x-x0)+y0, 'k')
# 
# 
# im.set_clim((-4.5,4.5))
# ax.set_ylabel('y [m]')
# ax.set_xlabel('x [m]')
# cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.63)
# #cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.79)
# cbar.ax.set_ylabel('h [mm]', rotation=270, labelpad=20)
# ax.set_ylim(y[-1], y[0])
# plt.tight_layout()
# plt.show()
# # plt.savefig('/home/samantha/Desktop/fig1.png', dpi=200)


mask_line = np.zeros_like(hmur)
xindex = np.arange(len(x))
yindex = -m*(xindex - idx_x0) + 685


for i in range(mask_line.shape[1]):
    yidx = int(yindex[i])
    mask_line[yidx,i]=1

# plt.figure()
# plt.imshow(mask_line)
# plt.show()
# 
# plt.figure()
# plt.imshow(hmur, cmap=colormap)
# plt.show()



mult = hmur*mask_line
linea = np.zeros_like(y, dtype='float')
for i in range(hmur.shape[0]):
    try:
        idx_value  = np.argwhere(mult[i,:])[0][0]
        linea[i] = mult[i,idx_value]
    except:
        pass


xtilde = y/np.sin(26.81*np.pi/180)
xtilde = xtilde[:615]
linea = linea[:615]

# plt.figure()
# plt.plot(xtilde, linea)
# plt.grid()
# plt.show()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


# for i in range(3, 12):
#     plt.figure()
#     plt.plot(running_mean(linea,i), label=str(i))
#     plt.title(str(i))
#     #plt.legend()
#     plt.grid()
#     plt.show()


mean_line = running_mean(linea, 5)
xtilde_mean  = running_mean(xtilde, 5)


peaks, _ = find_peaks(mean_line, height=0, distance = 50)

x1 = xtilde_mean[peaks[0]]
x2 = xtilde_mean[peaks[1]]


fig, ax = plt.subplots(figsize=(8.27,2))
ax.plot(xtilde_mean, mean_line, 'royalblue')
#ax.plot(xtilde_mean[peaks], mean_line[peaks], 'r.')
yarrow = 1.2
#ax.annotate(text='', xy=(x1,yarrow), xytext=(x2,yarrow), arrowprops=dict(arrowstyle='<->', overhang=0.8))
ax.annotate('', (x1,yarrow), (x2, yarrow), arrowprops=dict(facecolor='k', arrowstyle='<|-|>', mutation_scale=10, linewidth=1))
greek_letterz=[chr(code) for code in range(945,970)]
ax.text(-0.140, 1.27, greek_letterz[10])
ax.set_xlabel(r'$y/\sin\theta$ [m]')
ax.set_ylabel('$h$ [mm]')
ax.grid(alpha=0.6)
ax.set_ylim(ymax=1.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/wavelength.png')





long_de_onda = np.mean(np.diff(xtilde_mean[peaks]))
std_long_de_onda = np.std(np.diff(xtilde_mean[peaks]))







