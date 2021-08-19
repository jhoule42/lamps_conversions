# ******************************************************************************
#  Traitement des cartes de contributions
#  Auteur: Julien-Pierre Houle
# ******************************************************************************

import os
import numpy as np
import pandas as pd
import MultiScaleData as MSD
import glob
import matplotlib.pyplot as plt


# Create circular masks
def set_circle(a,r,v):
    x = np.arange(len(a)) - len(a)/2
    X,Y = np.meshgrid(x,x)
    a[(X*X+Y*Y) <= r*r] = v


def sum_circle(a,r):
    x = np.arange(len(a)) - len(a)/2
    X,Y = np.meshgrid(x,x)
    return np.sum(a[(X*X+Y*Y) <= r*r])


MSD_map = { scen: np.dstack([ MSD.Open(m)[0] for m in sorted(glob.glob(f'results/{scen}/*.hdf5')) ]) \
            for scen in os.listdir('results') if 'contribution_maps' not in scen}


scoto = np.loadtxt('public/Lights/scotopic.dat', skiprows=1)
wl, spectre = scoto.T
interval = (wl>= 400) & (wl < 600)

# Split sur 5 parties et faire la moyenne
moy_wl = [np.mean(x) for x in np.array_split(wl[interval], 5)]
moy_spectre = [np.mean(x) for x in np.array_split(spectre[interval], 5)]

scoto_maps = { key : np.dot(MSD_map[key], moy_spectre) * (moy_wl[1]-moy_wl[0]) for key in MSD_map }

# Save contribution maps for data manipulation
for keys in scoto_maps:
    np.save(f"scripts/np_contrib_maps/{keys}", scoto_maps[keys])


scens = {
    "fpub": "public",
    "fpubc": "public_conv",
    "fpri": "private",
    "fpric": "private_conv"}

for l in scens:
    contrib = np.load(f"scripts/np_contrib_maps/{scens[l]}.npy")
    x = [ sum_circle(contrib,r) for r in np.linspace(0,150,1001) ]
    plt.plot(np.linspace(0,150,1001),x,label=scens[l])

# Read conversions scenarios files
for conv in glob.glob("scenarios/scen_*.csv"):
    tab = pd.read_csv(conv)
    name = conv[10:-4]
    skys = dict()

    for l in scens:
        contrib = np.load(f"scripts/np_contrib_maps/{scens[l]}.npy")
        weight = np.zeros_like(contrib)

        for r,w in tab[["R",l]].to_numpy()[::-1]:
            set_circle(weight,r*10,w)

        # Intensite en unite naturel (3.75e-7: zenith naturel 5 deg angle vue)
        skys[l] = np.sum(weight*contrib) / 3.75e-7

    print(f"{name}:\t{sum(list(skys.values()))}")

# plt.imshow((np.log(scoto_map['private_conv'])))
