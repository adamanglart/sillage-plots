import numpy as np
from matplotlib import pyplot as plt
import scipy.io as sio
from figure_styling import *
plt.close('all')
plt.ion()


path = '/media/samantha/My Passport/these_sillage_plots/mediciones_wake/'
mediciones_folder_1 = 'mediciones-2023-06-22/'



data = sio.loadmat(path+mediciones_folder_1+'reflection.mat')

Fr       = data['Fr'][0]
Rs_mur   = data['Rs_mur'][0]
Rs_gra   = data['Rs_gra'][0]
eRs_mur  = data['eRs_mur'][0]
eRs_gra  = data['eRs_gra'][0]




fig, ax = plt.subplots(figsize=(7,3.5))
ax.plot(Fr , Rs_mur , '.--', c='royalblue', label='mur')
ax.errorbar(Fr , Rs_mur , yerr=eRs_mur , linestyle ='none', c='royalblue', capsize=1.7)
ax.plot(Fr , Rs_gra , '.-', c='firebrick', label='gra')
ax.errorbar(Fr , Rs_gra , yerr=eRs_gra , linestyle ='none', c='firebrick', capsize=1.7)
#ax.legend()
ax.text(0.77, 0.265, 'solid wall', color='royalblue')
ax.text(0.776, 0.1428, 'grating',    color='firebrick')
#ax.set_ylim([0,1])
ax.set_xlabel('Fr')
ax.set_ylabel('$|R|$')
ax.grid(alpha=0.6)
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/reflection_wake.eps', dpi=100)

# stop
# fig, ax = plt.subplots(figsize=(7,5))
# ax.plot(Fr , Rs_mur**2 , '.--', c='royalblue', label='mur')
# #ax.errorbar(Fr , Rs_mur , yerr=eRs_mur , linestyle ='none', c='royalblue', capsize=1.7)
# ax.plot(Fr , Rs_gra**2, '.-', c='firebrick', label='gra')
# #ax.errorbar(Fr , Rs_gra , yerr=eRs_gra , linestyle ='none', c='firebrick', capsize=1.7)
# #ax.legend()
# # ax.text(0.77, 0.265, 'solid wall', color='royalblue')
# # ax.text(0.776, 0.121, 'grating',    color='firebrick')
# ax.set_yscale('log')
# #ax.set_ylim([0,1])
# ax.set_xlabel('Fr')
# ax.set_ylabel('$|R|^2$')
# ax.grid(alpha=0.6)
# plt.show()
# 
# fig, ax = plt.subplots(figsize=(7,5))
# ax.plot(Fr , Rs_mur - Rs_gra , '.-', c='royalblue', label='mur')
# # ax.errorbar(Fr , Rs_mur , yerr=eRs_mur , linestyle ='none', c='royalblue', capsize=1.7)
# # ax.plot(Fr , Rs_gra , '.-', c='firebrick', label='gra')
# # ax.errorbar(Fr , Rs_gra , yerr=eRs_gra , linestyle ='none', c='firebrick', capsize=1.7)
# # #ax.legend()
# # ax.text(0.77, 0.265, 'solid wall', color='royalblue')
# # ax.text(0.776, 0.121, 'grating',    color='firebrick')
# # ax.set_ylim([0,1])
# ax.set_xlabel('Fr')
# ax.set_ylabel(r'$|\Delta R|$')
# ax.grid(alpha=0.6)
# plt.show()














