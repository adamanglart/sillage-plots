import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import quad, romberg, quadrature, fixed_quad, simpson
from numba import jit
from figure_styling import *
plt.close("all")
plt.ion()

# ---- matplotlib styling
fs = 13
rc('legend', fontsize=fs)
rc('axes',   labelsize=fs)
rc('axes',   labelsize=fs, titlesize=fs)
rc('xtick',  labelsize=fs)
rc('ytick',  labelsize=fs)
rc('lines',  markersize=7)
mpl.rcParams['axes.titlesize'] = fs



# Fr = 0.80
Frs = np.arange(0.4,1.5, 0.2)
CL = [300,200,100,50,50,50]

# fig, axs = plt.subplots(2,3,sharex=True, sharey=True, figsize=(14,10), dpi=100)
fig, axs = plt.subplots(2,3, figsize=(14,8))
ax = axs.ravel()


for kk in range(len(Frs)):
    Fr = Frs[kk]

    Lambda = 2 * np.pi * Fr**2
    
    # build mesh
    N = 200
    Nint = 5000
    
    # xmin = 0
    # xmax = 10
    xmin = -10
    xmax = 0

    
    ymin = 0
    ymax = 5
    
    Xvec = np.linspace(xmin, xmax, N) * Lambda
    Yvec = np.linspace(ymin, ymax, N) * Lambda
    Xg, Yg = np.meshgrid(Xvec, Yvec)
    
    # pressure distribution and integrand function
    
    
    @jit(nopython=False)
    def K0(theta):
        return 1 / (Fr**2 * np.cos(theta) ** 2)
    
    
    @jit(nopython=False)
    def P_hat(K):
        """
        Fourier transform of the pressure distribution.
        In this case, we are considering a gaussian pressure profile.
        """
        return np.exp(-(K**2) / (4 * np.pi) ** 2)
    
    
    @jit(nopython=False)
    def integrand(theta, X, Y):
        num = P_hat(K0(theta)) * np.exp(
            -1j * (-X * np.cos(theta) - Y * np.sin(theta)) / (Fr * np.cos(theta)) ** 2
        )
        den = Fr**4 * np.cos(theta) ** 4
        return num / den
    
    
    # calculate the integral
    Zfield = np.zeros((N, N), dtype="complex")
    
    for i in tqdm(range(N)):
        for j in range(N):
            # opcion 1, usar un integrador con funcion
            # Ztemp, err = quad(integrand, -np.pi, np.pi, args=(Xvec[i], Yvec[j]))
            # opcion 2, evaluar y usar trapecios
            theta_vals = np.linspace(-np.pi, np.pi, Nint, endpoint=True)
            vals = integrand(theta_vals, Xvec[i], Yvec[j])
            Ztemp = simpson(vals, theta_vals, dx=np.mean(np.diff(theta_vals)))
            Zfield[j, i] = 1j * np.pi * Ztemp
    
    Z = np.vstack((np.flipud(Zfield), Zfield))
    Z = np.imag(Z)
    
    xx = Xvec / Lambda
    yy = np.linspace(-ymax, ymax, 2 * N) / Lambda
    
    im = ax[kk].pcolormesh(xx,yy, Z, cmap = colormap)
    #plt.pcolormesh(Xvec / Lambda, np.linspace(-ymax, ymax, 2 * N) / Lambda, Z, cmap = colormap)
    #if kk==0 or kk==3:
    ax[kk].set_ylabel(r'$\tilde{Y}$')
    #if kk>2:
    ax[kk].set_xlabel(r'$\tilde{X}$')
    ax[kk].set_title('Fr = '+str("%.1f" % Fr))
    #plt.clim(-50,50)
    im.set_clim(-CL[kk], CL[kk])
    plt.colorbar(im, ax = ax[kk])
    #plt.tight_layout()
    #plt.show()
    
    #plt.savefig('/home/samantha/Desktop/sillage_Fr080.png')

plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/wakes_Fr.png')
#plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/wakes_Fr.eps', dpi=100)
