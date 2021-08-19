# ******************************************************************************
#  Graphiques des cartes de contributions
#  Auteur: Julien-Pierre Houle
# ******************************************************************************

import copy
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def set_circle(a,r,v):
    x = np.arange(len(a)) - len(a)/2
    X,Y = np.meshgrid(x,x)
    a[(X*X+Y*Y) <= r*r] = v

# Load contributions maps
dict_cartes = { m[:-4]: np.load(f"scripts/np_contrib_maps/{m}") for m \
                            in (os.listdir('scripts/np_contrib_maps'))}

dict_cartes.update({ 'total': (dict_cartes['public'] + dict_cartes['private']) })
dict_cartes.update({ 'total_conv': (dict_cartes['public_conv'] + dict_cartes['private_conv']) })

cb_labels = {'public':       'Contribution des sources de la bande scotopic [W/m^2/sr]',
             'public_conv':  "Contrib",
             'private':      "Contrib",
             'private_conv': "Contrib",
             'total': 'Contributions des sources [W/m^2/sr]',
             'total_conv': "Contrib" }

# Set nan values to black
cmap = copy.copy(plt.cm.get_cmap("inferno"))
cmap.set_bad(color='black')

for key in dict_cartes:
    set_circle(dict_cartes[key], 5, 0) # Set 0 to 500 meters radius
    plt.imshow(dict_cartes[key], cmap=cmap, norm=LogNorm(vmin=1e-11, vmax=1e-8))
    plt.colorbar(label=cb_labels[key])
    plt.axis('off')
    plt.savefig(f'results/contribution_maps/{key}', bbox_inches='tight')
    plt.close()



# plt.imshow(np.log(carte), cmap=cmap)
# plt.imshow(carte, cmap=cmap, norm=LogNorm())
# plt.colorbar(label='Contribution des sources de la bande scotopic [W/m^2/sr]')
# plt.axis('off')
# plt.savefig(f'results/contribution_maps/{m[:-4]}', bbox_inches='tight')
# plt.close()
