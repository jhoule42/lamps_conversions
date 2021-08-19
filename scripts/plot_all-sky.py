# Script d'analyse des resulats d'Illumina pour cartes du ciel
# Auteur : Julien-Pierre Houle

# Importer les librairies
import numpy as np
import matplotlib.pyplot as plt
import pytools as pt
import MultiScaleData as MSD
import hdftools as hdf
import os, sys
from math import log10, sqrt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
sys.path.append("git/Illumina-scripts/Sources")
from custum_cmap import standart_cmap


# Elements de la Matrice
m = np.zeros((6, 12, 10))  # Matrice de (el, az, wl)
elevation = [15, 30, 45, 60, 75, 90]
azimuth = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
wavelength = [397.0, 432.0, 467.0, 502.0, 537.0, 572.0, 607.0, 642.0, 677.0, 712.0]

listes = [elevation, azimuth, wavelength]
carte_sky_dict = {}
m_dict = {}


dict_SBbg = {"U": 22.03, "B": 22.73, "V": 21.93, "R": 21.18, "I": 20.03, "SQM": 22.0}  # (mag/arcsec^2)
dict_Rbg = {"U": 1.72e-07, "B": 2.05e-07, "V": 2.22e-07, "R": 5.10e-07, "I": 7.65e-07} # (W*sr^-1*m^-2)
dict_R0 = {"U": 111.8, "B": 254.3, "V": 131.4, "R": 151.2, "I": 78.7, "SQM": 270.0}    # (W*sr^-1*m^-2)

# SCÉNARIOS de conversions
path = ["git/Illumina-scripts/Results_3/Results/"]

light_source = ["HS", "private", "university"]
saisons = ["_ete", "_hiver"]
light_type = ["_2200K", "_3000K", "_actual", "_AMBR"]
txt = ["/data.txt"]
scenarios = []
nom_scenarios = []
saison_type = []
dict_sources_comb = {}
dict_comb_total = {}
R_dict = {}
Sky_comb_dict = {}


# Creer le path de chaque scenario
for p in path:
    for source in light_source:
        for saison in saisons:
            for type in light_type:
                for t in txt:
                    scenarios.append(p + source + saison + type + t)


# Creer le nom de chaque scenario
for source in light_source:
    for saison in saisons:
        for type in light_type:
            nom_scenarios.append(source + saison + type)


# Creer les sous-noms saison-light_type
for saison in saisons:
    for type in light_type:
        saison_type.append(saison + type)


for s_t in saison_type:
    scen_comb = []
    for scen in nom_scenarios:  # Pour tout les scenarios

        if s_t in scen:
            scen_comb.append(scen)
            dict_sources_comb.update({s_t: scen_comb})   # Dict saison-light_type : liste 3 light_source



# SPECTRES d'Integration
path_spct = "git/Illumina-scripts/Spectres"
scoto = path_spct + "/scotopic.dat"
sqm = path_spct + "/sqm.dat"
liste_spectres = [sqm]
nom_spectres = ['SQM']


print("\n------------------------------------------------")
for index1, scenario in enumerate(scenarios):    # Itérations a travers chaque scénarios

    R_total = []

    with open(scenario) as f:
        for line in f:

            line = line.split('-', 2)
            elevation = int(line[0][16:])
            azimuth = int(line[1][14:])
            wavelength = int(line[2][11:14])
            valeur = (line[2][17:])
            #print("{}".format(line))

            # Remplir la matrice pour scens individuels
            if elevation != 0:
                m[listes[0].index(elevation), listes[1].index(azimuth), listes[2].index(wavelength)] = valeur
                m[-1, :, :] = m[-1, 0, :]  # Lorsqu'on est au zénith
                m_dict.update({nom_scenarios[index1]: m.copy()})    # Dict de Matrice


    elevation = [15, 30, 45, 60, 75, 90]  # Ouin je triche un peu my bad
    azimuth = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]


    # Intégration sur les spectres
    for index2, spct in enumerate(liste_spectres):
        spectre = np.loadtxt(spct, skiprows=1)
        wl, spectre = spectre.T
        interval = (wl >= 380) & (wl <= 730)

        # Split sur 10 parties et fait la moyenne
        moy_wl = [np.mean(x) for x in np.array_split(wl[interval], 10)]
        moy_spectre = [np.mean(x) for x in np.array_split(spectre[interval], 10)]

        # Afficher wl moyenner sur la courbe
        #plt.plot(wl[interval], spectre[interval])
        #plt.plot(moy_wl, moy_scoto,'ko')

        # Intégration sur les longueurs d'ondes
        print("\n",str(index1+1) + ". Integration de " + str(nom_scenarios[index1]) + " sur " + str(nom_spectres[index2]))
        sky = (moy_wl[1]-moy_wl[0]) * np.dot(m, moy_spectre)

        # Caster dans dict + somme
        carte_sky_dict.update({(nom_scenarios[index1]): sky.copy()})
        total_carte2 = np.sum([np.sum(layer) for layer in carte_sky_dict[nom_scenarios[index1]]])
        print("\tTotal Sky: ", total_carte2)


for idx, scen in enumerate(dict_sources_comb.values()): # Iterer a travers les listes
    R = np.zeros((3, 6, 12))

    for idx1, x in enumerate(scen):  # Pour chaque scenario de la liste
        for idx2, j in enumerate(carte_sky_dict[x]):    # Pour chaque array de 12
            for idx3, valeur in enumerate(j):    # Pour chaque element des 12

                R[idx1][idx2][idx3] = valeur
    R_dict.update({(idx): R})   # Dictionnaire de array R


# Iterer a travers le R_dict pour faire le somme
for liste_source in R_dict.keys():
    print(liste_source)

    Sky_comb = (np.sum(R_dict[liste_source], axis = 0))
    print(Sky_comb)
    Sky_comb_dict.update({(liste_source): Sky_comb})


# Plot dans la fonction pour la carte du ciel
for idx, sky_map in enumerate(Sky_comb_dict.values()):
    pt.plot_allsky(azimuth, elevation, sky_map,
                   interp="None",
                   cmap="inferno",
                   autogain=False,
                   clabel = "Sky brightness for %s band [$W/sr/m^2$]" % nom_spectres[index2],
                   showpts=False,
                   log=True,
                   vmin = 3e-6,
                   vmax = 5e-9)
    plt.savefig("/home/jhoule42/Documents/Resultats_Sherbrooke3/Carte_Ciel/Ra/%s_log_none.png" % saison_type[idx], bbox_inches="tight")
    plt.close()

# plt.imshow(a[0])
# plt.imshow(a[0], cmap="inferno")

# Print ZENITHS
print("\tZenith : ", sky[-1][0])

elevation = [15, 30, 45, 60, 75, 90]  # Ouin je triche un peu my bad
azimuth = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]

newcmap = standart_cmap()
cmap = ListedColormap(newcmap(np.linspace(0.0, 1.0, 256)))

#delta_mesure = [1.7, 2.28 ]


# Convertir SKY BRIGHTNESS
print("\n------------------------------------------------")
for index1, sky_map in enumerate(Sky_comb_dict.values()):    # Itérations a travers scénarios

    print("Conversion to SB: ", saison_type[index1])

    # Convertir la brillance en magnitude
    layer = Sky_comb_dict[index1]
    Rbg = dict_R0["SQM"] * (10 ** (-0.4*dict_SBbg["SQM"]))
    SB_sky = (-2.5*np.log10((layer + Rbg)/dict_R0["SQM"]))
    print("\tZenith : ", SB_sky[-1][0], '\n')

#    SB_sky -= delta_mesure[index1]

    # Plot dans la fonction pour la carte du ciel
    pt.plot_allsky(azimuth, elevation, SB_sky,
                   interp="None",
                   cmap = cmap,
                   autogain=False,
                   showpts=False,
                   vmin = 14.0,
                   vmax = 24.0)

    cbar = plt.colorbar(ticks=list(range(14,25)))
    cbar.set_label("Sky brightness for %s band [$mag/arcsec^2$]" % nom_spectres[index2])
    cbar.ax.invert_yaxis()
    cbar.ax.set_yticklabels([str(i) for i in range(14,25)])


    plt.savefig("/home/jhoule42/Documents/Resultats_Sherbrooke3/Carte_Ciel/SB/%s_SB.png" % saison_type[index1], bbox_inches="tight")
    plt.close()

    plt.ion()   # Idealement ajouter dans le script d'initialisation
    plt.show()

print("------------------------------------------------")
print("Execution complete!")
