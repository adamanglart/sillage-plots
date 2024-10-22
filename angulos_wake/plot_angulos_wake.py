import numpy as np
from matplotlib import pyplot as plt
import sys
sys.path.insert(0, '/media/samantha/My Passport/these_sillage_plots/')
from figure_styling import *
from matplotlib import cm
plt.close('all')
plt.ion()


alpha_total = np.load('/media/samantha/My Passport/these_sillage_plots/angulos_wake/alpha.npy')
theta_total = np.load('/media/samantha/My Passport/these_sillage_plots/angulos_wake/theta.npy')

# v = 900, variando H 
# /media/samantha/SP PHD U3/mediciones-2023-06-16
Fr_variando_h    = [0.73,  0.59,  0.51,  0.46,  0.42,  0.39,  0.36,  0.34]
theta_variando_h = [26.58, 32,    33.74, 33.66, 35.66, 36.33, 36.66, 33.66]
error_theta_varh = [1.49,  3.48,  3.22,  2.81,  3.91,  2.86,  3.57,  4.64]

# H = 0.02, variando v
# /media/samantha/SP PHD U3/mediciones-2023-06-06
Fr_variando_v    = [0.73,  0.81,  0.89]
theta_variando_v = [26.58, 23.55, 20.88  ]
error_theta_varv = [1.49,  3.75, 2.44]


# H = 0.02, variando v cerca de v = 900
# /media/samantha/My Passport/mediciones-2023-06-22
Fr_variando_v9    = [0.70,  0.71,  0.72,  0.73,  0.74,  0.75,  0.76,  0.77,  0.78,  0.79,  0.80]
theta_variando_v9 = [26.66, 26.66, 29.33, 26.81, 28.44, 25.33, 26.66, 25.77, 20.88, 24.43, 24.89]
error_theta_varv9 = [2.44,  2.00,  2.21,  5.67,  4.38, 2.38,   3.68, 4.92,   1.22,  2.22, 1.44]




Frs_total = np.arange(0.01,4,0.01)

#plt.figure(figsize=(7,5))
fig, ax = plt.subplots(figsize=(7,5))

ax.plot(Frs_total, theta_total, '-', c='k', label=r'$\theta$')
ax.plot(Fr_variando_h, theta_variando_h, '.', c='firebrick')
#ax.errorbar(Fr_variando_h, theta_variando_h, linestyle ='none', yerr = error_theta_varh, c='firebrick', capsize=1.7)
ax.plot(Fr_variando_v, theta_variando_v, '.', c='firebrick')
#ax.errorbar(Fr_variando_v, theta_variando_v, linestyle ='none', yerr = error_theta_varv, c='firebrick', capsize=1.7)
ax.plot(Fr_variando_v9, theta_variando_v9, '.', c='firebrick')
#ax.errorbar(Fr_variando_v9, theta_variando_v9, linestyle ='none', yerr = error_theta_varv9, c='firebrick', capsize=1.7)
ax.set_xlabel('Fr')
ax.set_ylabel(r'$\theta$ [degrees]')
ax.set_xlim(Frs_total[0], Frs_total[-1])

# inset axes....
x1, x2, y1, y2 = 0.3, 0.95, 10, 45  # subregion of the original image
axins = ax.inset_axes(
    [0.5, 0.08, 0.45, 0.61],
    xlim=(x1, x2), ylim=(y1, y2))#, xticklabels=[], yticklabels=[])
    #xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
axins.plot(Frs_total, theta_total, '-', c='k', label=r'$\theta$')
axins.plot(Fr_variando_h, theta_variando_h, '.', c='firebrick')
axins.errorbar(Fr_variando_h, theta_variando_h, linestyle ='none', yerr = error_theta_varh, c='firebrick', capsize=1.7)
axins.plot(Fr_variando_v, theta_variando_v, '.', c='firebrick')
axins.errorbar(Fr_variando_v, theta_variando_v, linestyle ='none', yerr = error_theta_varv, c='firebrick', capsize=1.7)
axins.plot(Fr_variando_v9, theta_variando_v9, '.', c='firebrick')
axins.errorbar(Fr_variando_v9, theta_variando_v9, linestyle ='none', yerr = error_theta_varv9, c='firebrick', capsize=1.7)
axins.grid(True, alpha=0.6)

ax.indicate_inset_zoom(axins, edgecolor="gray")

#plt.xlim([0, 1.25])
#plt.ylim([0, 40])
#plt.legend(edgecolor='None', facecolor='None')
#ax.grid(True, alpha=0.6)
plt.show()
