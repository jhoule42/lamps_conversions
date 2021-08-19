# Scripts d'extraction du modèle Gambons
# Date: Avril 2021

import os
import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
N = 1000

fine_els = np.linspace(0,np.pi/2,N//2+1)
fine_azs = np.linspace(0,2*np.pi,2*N+1)

az_c, el_c = np.meshgrid(fine_azs[:-1] + np.pi/(2*N), fine_els[:-1] + np.pi/(2*N))

avg_L = np.zeros((N//2,N*2))

for idx, file in enumerate(os.listdir('natural/GAMBONS'),1):
    print(idx)
    print('Opening file')
    az, al, L = np.loadtxt(f'natural/GAMBONS/{file}', delimiter=', ').T
    az %= 360
    az = np.tile(np.deg2rad([az-360,az,az+360]).flatten(),2)
    el = np.deg2rad(90-al)
    el = np.tile(el,3)
    el = np.array([el,-el]).flatten()
    L = np.tile(L,6)

    pts = np.array([az,el]).T

    print('Interpolating')
    fine_L = griddata(pts, L, (az_c, el_c))

    avg_L += (fine_L-avg_L) / idx

np.save('natural/avg_L',avg_L)

ax = plt.subplot(111,polar=True)
plt.pcolormesh(fine_azs,fine_els,avg_L,cmap='inferno')
plt.colorbar(label='Intensité lumineuse (cd)')
plt.ylim(0,np.pi/2)
ax.set_theta_zero_location('N')
ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
