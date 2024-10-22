import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

from matplotlib import rc
import matplotlib as mpl

# ---- general 
import cmocean
colormap = cmocean.cm.gray

# ---- matplotlib styling
fs = 12
rc('legend', fontsize=12)
rc('axes',   labelsize=fs)
rc('axes',   labelsize=fs, titlesize=fs)
rc('xtick',  labelsize=fs)
rc('ytick',  labelsize=fs)
rc('lines',  markersize=7)
mpl.rcParams['axes.titlesize'] = fs

mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif'] = ['Times New Roman']+ mpl.rcParams['font.serif']
mpl.rcParams["mathtext.fontset"] = "stix"

from matplotlib import cm
plt.close('all')
plt.ion()


color1 = 'black'
color2 = 'firebrick'

#path = '/media/samantha/My Passport/mediciones-2023-01-11/'
# path = '/media/samantha/SP PHD U3/mediciones-2023-06-02/'

data_calib = np.load('speed_vs_real_speed_calib.npy')


ss = data_calib[0]
vv = data_calib[1]

ss = ss[:32]
vv = vv[:32]

idx_lineal = 17

z = np.polyfit(ss[:idx_lineal], vv[:idx_lineal], 1)
p = np.poly1d(z)

formula1 = 'v$_{\mathregular{real}} = $'+str(round(z[0],2))+'v$_{\mathregular{arduino}}$'+str(round(z[1],2))
## COMPARACION CON CALCULO DE LA VELOCIDAD

microsteps_por_vuelta = 400
R = 23.38 # mm # radio de la ruedita
v_arduino_teorica = (2*np.pi/microsteps_por_vuelta)*ss*R

z2 = np.polyfit(v_arduino_teorica[:idx_lineal], vv[:idx_lineal], 1)
p2 = np.poly1d(z2)

formula2 = 'v$_{\mathregular{real}} = $'+str(round(z2[0],2))+'v$_{\mathregular{arduino}}^{\mathregular{theo}}$'+str(round(z2[1],2))


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8,3.9))

ax1.plot(ss, vv, '.-', color=color1)
ax1.plot(ss[:idx_lineal], p(ss[:idx_lineal]), '-', label = formula1, color=color2)
ax1.text(700,225, formula1, color = color2)#, backgroundcolor = 'w')
ax1.text(170,465, '(a)', color = 'k')#, backgroundcolor = 'w')

ax2.plot(v_arduino_teorica, vv, '.-', color=color1)
ax2.plot(v_arduino_teorica[:idx_lineal], p2(v_arduino_teorica[:idx_lineal]), '-', label = formula2, color=color2)
ax2.text(p(700),225, formula2, color = color2)#, backgroundcolor = 'w')
ax2.text(68,465, '(b)', color = 'k')#, backgroundcolor = 'w')

ax1.grid(alpha=0.6)
ax2.grid(alpha=0.6)

ax1.set_ylabel(r'v$_{\mathregular{real}}$ [mm/s]')
ax1.set_xlabel(r'v$_{\mathregular{arduino}}$ [microsteps/s]', )
ax2.set_xlabel(r'v$_{\mathregular{arduino}}^{\mathregular{theo}}$ [mm/s]', )

plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/fit_calibracion.eps')





## QUÉ NUMEROS DE FROUDE PODEMOS ALCANZAR

g = 9.82 # m/s^2

def froude(v,h):
    v = v/1000
    return v/np.sqrt(g*h)

hs = np.arange(10,100,10)

cmap  = cm.get_cmap('tab20b')
cmap  = cm.get_cmap('Blues')
markers = ['o', 'h','v', 's', 'D', 'P', '>', 'X', '<', 'd']
ms = 3


# plt.figure()
# for h in hs:
#     plt.plot(vv, froude(vv, h/1000), '.-', label= str(h), c = cmap(h/100))
# plt.xlabel(r'v$_{\mathregular{real}}$ [mm/s]')
# #plt.ylabel('Fr = '+r'$\frac{\mathregular{s}}{\sqrt{gh}}$')
# plt.ylabel('Fr')
# legend = plt.legend()
# legend.set_title("h [mm]")   
# plt.grid()
# plt.tight_layout()
# plt.show()

plt.figure(figsize=(8,4.5))

for i in range(len(hs)):
    h = hs[i]
    plt.plot(ss, froude(vv, h/1000),  marker = markers[i], markersize = ms, markeredgewidth=1.5, markerfacecolor='white', label= str(h),  c =  cmap((len(hs)-i+2)/(len(hs)+1)))

plt.xlabel(r'v$_{\mathregular{arduino}}$ [microsteps/s]', )
#plt.ylabel('Fr = '+r'$\frac{\mathregular{s}}{\sqrt{gh}}$')
plt.ylabel('Fr')
#legend = plt.legend()
legend = plt.legend(loc=2, bbox_to_anchor=(1.02,0.9))
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
legend.set_title("$h$ [mm]")   
plt.grid()
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/froude_calibracion.eps')





