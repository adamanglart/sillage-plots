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
from scipy.io import savemat
plt.close('all')
plt.ion()

def R_from_data(kx_calculated, ky_calculated, y1, y2, eta1, eta2, x_vertical):
    """
    Calculates the reflection coefficient from the field at two points: (x0, y1) and (x0, y2), assuming that we know kx and ky.
    """
    num = (eta2*np.exp(1j*ky_calculated*y1) - eta1*np.exp(1j*ky_calculated*y2))*np.exp(-1j*kx_calculated*x_vertical)
    den = np.exp(1j*ky_calculated*y1)*np.exp(-1j*ky_calculated*y2) - np.exp(-1j*ky_calculated*y1)*np.exp(1j*ky_calculated*y2)
    R_calculated = num/den
    return R_calculated



def obtain_reflection_coefficient(height, idx_vertical_0, idx_vertical_f, idx_first, idx_last, delta_idx):



    x = np.linspace(0,  height.shape[1]*pixel_size, height.shape[1])
    y = np.linspace(0, -height.shape[0]*pixel_size, height.shape[0])

    # I take several vertical lines and I look for ky

    vertical_lines = height[:,idx_vertical_0:idx_vertical_f]
    ft_vertical_lines = np.fft.fftshift(np.fft.fft(vertical_lines, axis=0), axes=0)
    ft_vertical_mean = np.mean(np.abs(ft_vertical_lines), axis=1)
    ky_fft  = np.fft.fftshift(np.fft.fftfreq(len(y), np.diff(y)[0]))*2*np.pi

    idx0 = len(ft_vertical_mean)//2 + 1
    idx_ky1 = np.argmax(ft_vertical_mean[idx0:]) + idx0
    ky_calculated = np.abs(ky_fft[idx_ky1])

    # plt.figure()
    # plt.plot(ky_fft, ft_vertical_mean)
    # plt.plot(ky_fft[-idx_ky], ft_vertical_mean[idx_ky], 'ko')
    # plt.xlim([-500,500])
    # plt.title('fft in y')
    # plt.grid()
    # plt.show()
    
    
    # plt.figure()
    # iis = np.arange(0, 600, 100)
    # for i in iis:
    #     plt.plot(ky_fft, np.abs(ft_vertical_lines[:,i]), label=str(i))
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    
    # Now I know ky, at least for that vertical line
    
    # I take several horizontal lines and I look for kx
    
    rectangle_lines  = height[idx_first:idx_last + delta_idx,:]
    ft_rect = np.fft.fftshift(np.fft.fft(rectangle_lines, axis=1), axes=1)
    ft_mean = np.mean(np.abs(ft_rect), axis=0)
    kx_fft  = np.fft.fftshift(np.fft.fftfreq(len(x), np.diff(x)[0]))*2*np.pi
    idx0 = len(ft_mean)//2 + 3
    idx_kx1 = np.argmax(ft_mean[idx0:]) + idx0
    kx_calculated = kx_fft[idx_kx1]
    
    # plt.figure()
    # plt.plot(kx_fft, ft_mean)
    # plt.plot(kx_fft[idx_kx1], ft_mean[idx_kx1], 'ko')
    # plt.xlim([-500,500])
    # plt.title('Rectangle\nfft in x')
    # plt.grid()
    # plt.show()
    
    
    vertical_idxs = np.arange(idx_vertical_0, idx_vertical_f)
    idxs1 = np.arange(idx_first,idx_last)
    R  = np.zeros((len(vertical_idxs), len(idxs1)),  dtype=complex)
    for ii in range(len(vertical_idxs)):
        idx_vertical = vertical_idxs[ii]
        x_vertical   = x[idx_vertical]
        for jj in range(len(idxs1)):
            idx1 = idxs1[jj]
            idx2 = idx1 + delta_idx
        
            y1 = y[idx1]
            y2 = y[idx2]
        
            eta1 = height[idx1, idx_vertical]
            eta2 = height[idx2, idx_vertical]
        
            R_calculated = R_from_data(kx_calculated, ky_calculated, y1, y2, eta1, eta2, x_vertical)
        
            R[ii,jj] = R_calculated

    return R


path = '/media/samantha/My Passport/these_sillage_plots/mediciones_wake/'
mediciones_folder = 'mediciones-2023-06-22/'
vs = np.arange(870,1010,10)





g = 9.82    # m/s^2
h = 0.02    # m


velocidad_dim = 0.36*vs - 2.55 # mm/s
Fr            = (velocidad_dim/1000)/np.sqrt(g*h)

Rs_mur = np.zeros_like(vs, dtype=float)
Rs_gra = np.zeros_like(vs, dtype=float)
eRs_mur = np.zeros_like(vs, dtype=float)
eRs_gra = np.zeros_like(vs, dtype=float)


for ii in range(len(vs)):


    height_mur = np.load(path + mediciones_folder + ' mur_'+str(vs[ii])+'.npy')
    height_gra = np.load(path + mediciones_folder + ' gra_'+str(vs[ii])+'.npy')
    
    
    velocidad_dim = 0.36*vs[ii] - 2.55 # mm/s
    froude_h      = (velocidad_dim/1000)/np.sqrt(g*h)
    
    print(r'Fr = %.2f' % froude_h)
    print('v = %.2f' % (velocidad_dim/10) +' cm/s')
    print('v = %d' % vs[ii] +' ms/s')
    
    ftp_proc_parameters = input_output.read_parameter_file(path+'processing_parameters.yaml')
    pixel_size  = ftp_proc_parameters['MEASUREMENT']['pixel_size']
    
    hmur = np.flipud(height_mur)
    hgra = np.flipud(height_gra)
    hmur  = hmur[:,:1750]
    hgra  = hgra[:,:1750]
    
    # Vertical
   #  idx_vertical_0 = 600
   #  idx_vertical_f = 1000
    idx_vertical_0 = 600
    idx_vertical_f = 640
    # horizontal
    idx_first = 100
    #idx_last  = 200
    idx_last  = 300
    #delta_idx = 52
    delta_idx = 70
    
    vertical_idxs = np.arange(idx_vertical_0, idx_vertical_f)
    idxs1 = np.arange(idx_first,idx_last)
    
    Rmur = obtain_reflection_coefficient(hmur, idx_vertical_0, idx_vertical_f, idx_first, idx_last, delta_idx)
    Rgra = obtain_reflection_coefficient(hgra, idx_vertical_0, idx_vertical_f, idx_first, idx_last, delta_idx)
    
    # plt.figure()
    # plt.pcolor(idxs1, vertical_idxs, abs(Rmur), cmap=cmocean.cm.ice)
    # plt.title('mur')
    # plt.colorbar()
    # plt.grid()
    # plt.show()
    # 
    # plt.figure()
    # plt.pcolor(idxs1, vertical_idxs, abs(Rgra), cmap=cmocean.cm.ice)
    # plt.title('gra')
    # plt.colorbar()
    # plt.grid()
    # plt.show()
    
    Rmean_mur = np.mean(abs(Rmur), axis=1)
    Rmean_gra = np.mean(abs(Rgra), axis=1)
    
    mR_mur = np.mean(Rmean_mur)
    eR_mur = np.std(Rmean_mur)
    
    mR_gra = np.mean(Rmean_gra)
    eR_gra = np.std(Rmean_gra)
    
    mur = "Rmur = {0:.2f} +- {1:.2f} ".format(mR_mur, eR_mur) 
    gra = "Rgra = {0:.2f} +- {1:.2f} ".format(mR_gra, eR_gra) 
    # print(mur)
    # print(gra)
    
    plt.figure(figsize=(6, 3))
    plt.plot(vertical_idxs, Rmean_mur, label=mur)
    plt.plot(vertical_idxs, Rmean_gra, label=gra)
    # plt.text(600, 0.7,  gra)
    # plt.text(600, 0.65,  mur)
    plt.legend()
    plt.ylim([0,1])
    plt.ylabel('$|R|$')
    plt.xlabel('idx y')
    plt.title(vs[ii])
    plt.grid()
    plt.tight_layout()
    plt.show()
    
    Rs_mur[ii]  = mR_mur
    Rs_gra[ii]  = mR_gra
    eRs_mur[ii] = eR_mur
    eRs_gra[ii] = eR_gra
    

plt.figure()
plt.plot(Fr, Rs_mur, '.-', label='mur')
plt.plot(Fr, Rs_gra, '.-', label='gra')
plt.errorbar(Fr, Rs_mur, yerr=eRs_mur,  label='mur')
plt.errorbar(Fr, Rs_gra, yerr=eRs_gra,  label='gra')
plt.ylim([0, 1])
plt.legend()
plt.grid()
plt.show()
plt.savefig(path+mediciones_folder+'reflection.png')


mdic = {"Fr": Fr, "Rs_mur": Rs_mur, "eRs_mur": eRs_mur, "Rs_gra": Rs_gra, "eRs_gra": eRs_gra}
savemat(path+mediciones_folder+'reflection.mat', mdic)

