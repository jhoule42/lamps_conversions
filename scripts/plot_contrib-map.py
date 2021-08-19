# Script d'analyse des resulats d'Illumina pour cartes de contribution
# Auteur : Julien-Pierre Houle
# Note : Executer a parit de HOME. Les changement de dossier peuvent etre modifie au besoin.

# Importer les librairies
import numpy as np
import matplotlib.pyplot as plt
import pytools as pt
import MultiScaleData as MSD
import hdftools as hdf
from glob import glob
import os
from math import log10, sqrt

# Parametres
r_max = 7500*sqrt(2)   # rayon max entre observateur-limite domaine (m)
scen_comparaison = ["Actu_ete", "Actu_hiver"]
scen_ete = ["2200k_ete", "2700k_ete", "3000k_ete", "Actu_ete", "Ambree_ete"]
scen_hiver = ["2200k_hiver", "2700k_hiver", "3000k_hiver", "Actu_hiver", "Ambree_hiver"]

nom_scenarios  = ["2200k_ete", "2200k_hiver", "2700k_ete", "2700k_hiver", "3000k_ete",
                  "3000k_hiver", "Actu_ete", "Actu_hiver", "Ambree_ete", "Ambree_hiver"]


# Initialiser un objet MSD batard vide
path1 = "git/Illumina-scripts/Results_Sherb_R2"
path2 = "git/Illumina-scripts/Results_3/Results/"
carte_total = MSD.Open(path1 + "/Actu_ete/elevation_angle_90-azimuth_angle_0-wavelength_397.0.hdf5")
carte_total_diff = MSD.Open(path1 + "/Actu_ete/elevation_angle_90-azimuth_angle_0-wavelength_397.0.hdf5")


light_source = ["HS", "private", "university"]
saisons = ["_ete", "_hiver"]
light_type = ["_2200K", "_3000K", "_actual", "_AMBR"]
scenarios = []
saison_type = []


for source in light_source:
    for saison in saisons:
        for type in light_type:
                scenarios.append(source + saison + type)

# Creer les sous-noms saison-light_type
for saison in saisons:
    for type in light_type:
        saison_type.append(saison + type)


carte_dict = {}
carte_dict_spct = {}
os.chdir(path2)


# PLOT CARTES CONTRIBUTIONS
# Iteration a travers chaque scenario
print("\n-----------------------------------------")
for index, scen in enumerate(scenarios):
    os.chdir(scenarios[index])

    # Obtenir dict longueurs d'ondes : objet MSD
    hdf5_files = sorted(glob("*.hdf5"))
    MSD_dict = {int(s.split('_')[-1].split('.')[0]): MSD.Open(s) for s in hdf5_files}
    print("MSD directory created: " + scenarios[index])
    os.chdir("../..")

    # SPECTRES d'Integration
    path_spct = "Spectres"
    scoto = path_spct + "/scotopic.dat"
    JC_U = path_spct + "/JC_U.dat"
    JC_B = path_spct + "/JC_B.dat"
    JC_V = path_spct + "/JC_V.dat"
    JC_R = path_spct + "/JC_R.dat"
    JC_I = path_spct + "/JC_I.dat"
    sqm = path_spct + "/sqm.dat"

    liste_spectres = [sqm]
    nom_spectres = ["SQM"]
    os.chdir("..")

    # Iteration a travers chaque spectre
    for index2, spct in enumerate(liste_spectres):

        # Reset carte_total a chaque debut iteration
        for i, layer in enumerate(carte_total):
            carte_total[i][:] = 0

        spectre = np.loadtxt(spct, skiprows=1)
        wl, spectre = spectre.T
        interval = (wl>= 380) & (wl <= 730)

        # Split sur 10 parties et faire la moyenne
        moy_wl = [np.mean(x) for x in np.array_split(wl[interval], 10)]
        moy_spectre = [np.mean(x) for x in np.array_split(spectre[interval], 10)]

        for wl in sorted(MSD_dict.keys()):
            for i, layer in enumerate(MSD_dict[wl]):
                carte_total[i] += layer * moy_spectre[i] * (moy_wl[1]-moy_wl[0])

        # Creer un dict (scen, spct) : carte_contrib.
        carte_dict.update({(scen): carte_total.copy()})

        # Graphique contribution
        hdf.plot(carte_total, area = True, cmap="inferno", vmax = 5e-8)
        plt.axis("off")
        #plt.xlim(-8, 8)
        #plt.ylim(-8, 8)
        plt.colorbar(label = "Pixel contribution for %s band [$W/sr/m^2$]" % nom_spectres[index2])
        plt.savefig("/home/jhoule42/Documents/Resultats_Sherbrooke3/Carte_contrib/%s_%s.png" % (scenarios[index], nom_spectres[index2]), bbox_inches="tight")
        plt.close()

    os.chdir("Results_3/Results")
print("Contributions maps saved\n")



# Plot  diff CARTES CONTRIBUTIONS
# print("\n------------------------------------")
# carte_diff_dict = {}
# for comparaison in scen_comparaison:
#     for scen in sorted(carte_dict.keys()):
#
#         # Ne pas soustraire 2 cartes identiques
#         if scen not in scen_comparaison and comparaison.split("_")[1] in scen:
#             for i, layer in enumerate(carte_dict[scen]):
#
#                 # Calcul de la difference
#                 carte_total_diff[i] = carte_dict[comparaison][i][:] - carte_dict[scen][i][:]
#             carte_diff_dict[scen] = carte_total_diff.copy()
#
#             # Graphique diff contribution
#             hdf.plot(carte_total_diff, cmap="inferno", area=True, vmax = 70e-8)
#             plt.axis("off")
#             plt.colorbar(label = "Contribution (W/sr/m²)")
#             plt.savefig("/home/jhoule42/Documents/Resultats_Sherb_R2/Diff_Carte_contrib/diff_(%s-%s).png" % (comparaison, scen), bbox_inches="tight")
#             plt.close()
#             print("Diff contib map created: " + comparaison + " - " + scen)
# print("Diff Contributions maps saved")



# Integration sur le RAYON
def sum_circle(ds, radii):
    total = np.zeros_like(radii)
    for i in range(len(ds)):
        ny,nx = ds[i].shape
        Y0, X0 = (ny-1)/2, (nx-1)/2
        R = radii / ds.pixel_size(i)
        Y, X = np.ogrid[:ny,:nx]
        d2 = (X-X0)**2 + (Y-Y0)**2
        total += np.sum(ds[i][d2 <= R**2]) # A verifier
    return total


# Listes des differents rayons [start, stop, nb]
rayon = np.linspace(0, r_max, 251)

scen_ete = ["2200k_ete", "2700k_ete", "3000k_ete", "Ambree_ete"]
scen_hiver = ["2200k_hiver", "2700k_hiver", "3000k_hiver", "Ambree_hiver"]

# Somme cartes total
total_actu_ete = np.sum([np.sum(layer) for layer in carte_dict["Actu_ete"]])
total_actu_hiver = np.sum([np.sum(layer) for layer in carte_dict["Actu_hiver"]])
total_actu = [total_actu_ete, total_actu_hiver]


# GRAPHIQUE CONTRIBUTION RAYON (%)
print("\n--------------------------------------")
print("Graphiques Contribution (%):")
couleur = ['orange', 'green', 'blue', 'red']

for index3, saison in enumerate(scen_comparaison):
    print("\n Total ",saison)
    color = 0

    # Enlever pour graph contrib actu
    for index4, scen in enumerate(carte_dict.keys()):
        # Enlever pour graph contib actu
        if scen not in scen_comparaison and saison.split("_")[1] in scen:
            print(scen)
            sum_r_conv = []
            for r , radii in enumerate(rayon):
                sum_r_conv.append(sum_circle(carte_diff_dict[scen], radii))
                #sum_r_conv.append(sum_circle(carte_dict[saison], radii))


            # Convertir la reduction en pourcentage
            conv_pourc = [elem/total_actu[index3]*100 for elem in sum_r_conv]

            # Graphique de la reduction de brillance (%) selon le rayon
            plt.plot(rayon, conv_pourc, color = couleur[color])
            plt.xscale("log")
            plt.xlim(400, 11000)
            plt.xlabel("Radius (m)")
            plt.ylabel("Relative reduction of the artificial SQM radiance (%)")
            plt.legend(("2200k", "2700k", "3000k", "1800k"), loc="upper left")
            #plt.legend(("Summer", "Winter"), loc="upper left")
            color += 1
    plt.savefig("/home/jhoule42/Documents/Resultats_Sherb_R2/Graph_rayon/reduction_%s.png" %(saison), bbox_inches="tight")
    plt.close()


# GRAPHIQUE INTENSITY REDUCTION (W*sr^-1*m^-2)
print("\n------------------------------------")
print("Graphiques Intensity Reduction: \n")
for scen in sorted(carte_diff_dict.keys()):
    print(scen)
    result = []

    for r , radii in enumerate(rayon):
        result.append(sum_circle(carte_diff_dict[scen], radii))

    plt.plot(rayon, result)
    plt.xlabel("Radius (m)")
    plt.ylabel("Reduction intensity (W/m²/sr)")
    plt.legend(("2200k ete", "2200k hiver", "2700k ete", "2700k hiver", "3000k ete", "3000k hiver", "Amber ete", "Amber hiver"), loc="upper left")
plt.savefig("/home/jhoule42/Documents/Resultats_Sherb_R2/Graph_rayon/intensity_reduction.png", bbox_inches="tight")
plt.close()




# GRAPHIQUE SB SELON RAYON
print("\n------------------------------------")
print("Conversion to SB: \n")

# Somme cartes total
total_actu_ete = np.sum([np.sum(layer) for layer in carte_dict["Actu_ete"]])
total_actu_hiver = np.sum([np.sum(layer) for layer in carte_dict["Actu_hiver"]])
total_actu = [total_actu_ete, total_actu_hiver]


dict_SBbg = {"U": 22.03, "B": 22.73, "V": 21.93, "R": 21.18, "I": 20.03, "SQM": 22.0}  # (mag/arcsec^2)
dict_Rbg = {"U": 1.72e-07, "B": 2.05e-07, "V": 2.22e-07, "R": 5.10e-07, "I": 7.65e-07} # (W*sr^-1*m^-2)
dict_R0 = {"U": 111.8, "B": 254.3, "V": 131.4, "R": 151.2, "I": 78.7, "SQM": 270.0}    # (W*sr^-1*m^-2)

# ------------------------------
# Valeur modelisation (r = 0):  |
# (ete = 21.28 / hiver = 20.35) |
# ------------------------------

delta_mesure = [1.43, 2.05] # (ete / hiver)
couleur = ['orange', 'green', 'blue', 'red']

for index4, saison in enumerate(scen_comparaison):
    color = 0
    for k, scen in enumerate(sorted(carte_dict.keys())):
        SB = []
        if scen not in scen_comparaison and saison.split("_")[1] in scen:
            print(scen)

            # Incrementer de 1 l'indice de la couleur dans la liste

            for r , radii in enumerate(rayon):

                cercle_soustraire = (sum_circle(carte_dict[saison], radii))
                Radiance_a = (sum_circle(carte_dict[scen], radii)) # Somme du cercle conv

                exterior_circle = (total_actu[index4] - cercle_soustraire)
                total_conv = (Radiance_a + exterior_circle)

                # Convertir la brillance en magnitude
                Rbg = dict_R0["SQM"] * (10 ** (-0.4*dict_SBbg["SQM"]))
                SB.append(-2.5*log10((total_conv + Rbg)/dict_R0["SQM"]))

            for x, valeur in enumerate(SB):
                SB[x] -= delta_mesure[index4]

            plt.plot(rayon, SB, color = couleur[color])
            plt.xlabel("Radius (m)")
            plt.ylabel("SB (mag/arcsec^2)")
            plt.legend(("2200k", "2700k", "3000k", "1800k"), loc="upper left")
            color += 1
    plt.savefig("/home/jhoule42/Documents/Resultats_Sherb_R2/Graph_rayon/scen_magn_SQM_%s_modif.png" % (saison), bbox_inches="tight")
    plt.close()


# Retourner au path HOME
os.chdir("../../..")
print("\n------------------------------------")
print("Execution complete!")
