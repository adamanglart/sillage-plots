import numpy as np
import matplotlib.pyplot as plt
from figure_styling import *


plt.close("all")
plt.ion()

# All magnitudes in SI units

U = 0.4
g = 9.8

U = 1
g= 1

Npts = 1000
psi = np.linspace(0, np.pi / 2, Npts)
# psi = np.linspace(0, np.pi / 4, Npts)
# psi = np.linspace(np.pi / 4, np.pi / 2, Npts)

#plt.figure()
fig, ax = plt.subplots()
#for i in range(0, 10):
for i in range(0, 7):
    theta = i
    #x = (U**2) * theta / g * np.cos(psi) * (1 + np.sin(psi) ** 2)
    x = (U**2) * theta / g * np.cos(psi) * (1 -  0.5*np.cos(psi) ** 2)
    y = (U**2) * theta / g * np.cos(psi) ** 2 * 0.5*np.sin(psi)

    x= -x
    y = -y
    ax.plot(x, y, "k")
    ax.plot(x, -y, "k")

ax.axis("equal")
#plt.axis("off")
ax.grid(True, alpha = 0.6)
ax.set_xlabel(r"$\frac{xg}{U^2}$", fontsize=14)
ax.set_ylabel(r"$\frac{yg}{U^2}$", fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# cone
# xmax = np.ceil(np.max(x))
ax.plot(x, np.tan(19.5 * np.pi / 180) * x,  "--", color='firebrick')
ax.plot(x, -np.tan(19.5 * np.pi / 180) * x, "--", color='firebrick')
plt.tight_layout()
plt.show()
plt.savefig('/home/samantha/Dropbox/PhD/these/figures/wake/ship_wake.eps', dpi=100)
