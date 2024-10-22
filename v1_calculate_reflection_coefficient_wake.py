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


fig, ax = plt.subplots()
#im = ax.imshow(height,aspect='equal', cmap = colormap,  extent=[0, x[-1], y[-1], y[0]])
im = ax.imshow(height,aspect='equal', cmap = colormap)#,  extent=[0, x[-1], y[-1], y[0]])
im.set_clim((-4.5,4.5))
ax.set_ylabel('$y$ [m]')
ax.set_xlabel('$x$ [m]')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.63)
#cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink = 0.79)
cbar.ax.set_ylabel('h [mm]', rotation=270, labelpad=20)

plt.tight_layout()
plt.show()

def R_from_data(kx_calculated, ky_calculated, y1, y2, eta1, eta2):
    """
    Calculates the reflection coefficient from the field at two points: (x0, y1) and (x0, y2), assuming that we know kx and ky.
    """
    num = (eta2*np.exp(1j*ky_calculated*y1) - eta1*np.exp(1j*ky_calculated*y2))*np.exp(-1j*kx_calculated*x_vertical)
    den = np.exp(1j*ky_calculated*y1)*np.exp(-1j*ky_calculated*y2) - np.exp(-1j*ky_calculated*y1)*np.exp(1j*ky_calculated*y2)
    R_calculated = num/den
    return R_calculated


# I take a vertical line (fixed x) and I look for ky

idx_vertical = 1000
x_vertical   = x[idx_vertical]

vertical_line = height[:,idx_vertical]

ft_vertical_line = np.fft.fftshift(np.fft.fft(vertical_line))
ky_fft  = np.fft.fftshift(np.fft.fftfreq(len(y), np.diff(y)[0]))*2*np.pi
idx0 = len(ft_vertical_line)//2 + 3
idx_ky = np.argmax(np.abs(ft_vertical_line[idx0:])) + idx0
ky_calculated = np.abs(ky_fft[idx_ky])

plt.figure()
plt.plot(ky_fft, np.abs(ft_vertical_line))
plt.plot(ky_fft[-idx_ky], np.abs(ft_vertical_line[idx_ky]), 'ko')
plt.xlim([-500,500])
plt.title('fft in y')
plt.grid()
plt.show()


# Now I know ky, at least for that vertical line

# I take a fixed horizontal line (fixed y) and I look for kx


idx1 = 200
y1 = y[idx1]

line1 = height[idx1,:]

ft_line1 = np.fft.fftshift(np.fft.fft(line1))
kx_fft  = np.fft.fftshift(np.fft.fftfreq(len(x), np.diff(x)[0]))*2*np.pi
idx0 = len(ft_line1)//2 + 3
idx_kx1 = np.argmax(np.abs(ft_line1[idx0:])) + idx0
kx_calculated_1 = kx_fft[idx_kx1]

plt.figure()
plt.plot(kx_fft, np.abs(ft_line1))
plt.plot(kx_fft[idx_kx1], np.abs(ft_line1[idx_kx1]), 'ko')
plt.xlim([-500,500])
plt.title('Line1\nfft in x')
plt.grid()
plt.show()

# I choose anothe horizontal line

idx2 = 250
y2 = y[idx2]
line2 = height[idx2,:]

ft_line2 = np.fft.fftshift(np.fft.fft(line2))
kx_fft  = np.fft.fftshift(np.fft.fftfreq(len(x), np.diff(x)[0]))*2*np.pi
idx0 = len(ft_line2)//2 + 3
idx_kx2 = np.argmax(np.abs(ft_line2[idx0:])) + idx0
kx_calculated_2 = kx_fft[idx_kx2]

plt.figure()
plt.plot(kx_fft, np.abs(ft_line2))
plt.plot(kx_fft[idx_kx2], np.abs(ft_line2[idx_kx2]), 'ko')
plt.xlim([-500,500])
plt.title('Line2\nfft in x')
plt.grid()
plt.show()


print(ky_calculated)
print(kx_calculated_1)
print(kx_calculated_2)

kx_calculated = np.mean([kx_calculated_1, kx_calculated_2])


eta1 = height[idx1, idx_vertical]
eta2 = height[idx2, idx_vertical]


R_calculated = R_from_data(kx_calculated, ky_calculated, y1, y2, eta1, eta2)

print(np.abs(R_calculated))
