# Script pour plotter les
# Auteur : Julien-Pierre Houle

# Importer les librairies
import numpy as np
import matplotlib.pyplot as plt
import pytools as pt
import MultiScaleData as MSD
import hdftools as hdf

path = "git/Illumina-scripts/Spectres"

# SPECTRES d'Integration
scoto = path + "/scotopic.dat"
JC_U = path + "/JC_U.dat"
JC_B = path + "/JC_B.dat"
JC_V = path + "/JC_V.dat"
JC_R = path + "/JC_R.dat"
JC_I = path + "/JC_I.dat"
SQM = path + "/sqm.dat"

spct = [scoto, JC_U, JC_B, JC_V, JC_R, JC_I, SQM]
spct_list = [np.loadtxt(x, skiprows=1) for x in spct]

# Normaliser scoto sur 100
for wv in spct_list[0]:
    wv[1] *= 100

for wv in spct_list[-1]:
    wv[1] *= 100
    
# Plot le graphique
[plt.plot(i[:,0], i[:,1]) for i in spct_list]
plt.xlabel("wavelength (nm)")
plt.ylabel("Transmittance (%)")
plt.legend(("Scotopic", "U", "B", "V", "R", "I", "SQM"), loc="upper left")
plt.savefig("/home/jhoule42/Documents/Resultats_Sherbrooke/Spct/Spectres.png", bbox_inches="tight")
