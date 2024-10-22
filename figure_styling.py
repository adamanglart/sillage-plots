from matplotlib import rc
import matplotlib as mpl

# ---- general 
import cmocean
colormap = cmocean.cm.balance

# ---- matplotlib styling
fs = 12
rc('legend', fontsize=fs)
rc('axes',   labelsize=fs)
rc('axes',   labelsize=fs, titlesize=fs)
rc('xtick',  labelsize=fs)
rc('ytick',  labelsize=fs)
rc('lines',  markersize=7)
mpl.rcParams['axes.titlesize'] = fs

mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif'] = ['Times New Roman']+ mpl.rcParams['font.serif']
mpl.rcParams["mathtext.fontset"] = "stix"



