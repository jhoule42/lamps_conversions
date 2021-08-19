import MultiScaleData as MSD
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


sim_viirs = np.load('makeVIIRS/sim_viirs.npy')
sim_viirs = gaussian_filter(sim_viirs, 1, mode='constant')
viirs = MSD.Open('makeVIIRS/viirs.hdf5')
diff_viirs = viirs[0]-sim_viirs
sim_total = 1 - (sim_viirs/viirs[0])


# plt.figure()
# plt.imshow(sim_viirs, vmin=0, vmax=60)
# plt.colorbar(label = "Intensite lumineuse [$nW/cm^2/sr$]")# plt.title('Simulated')
# plt.axis('off')
# plt.savefig('makeVIIRS/Results/viirs_simulated', dpi=300, bbox_inches="tight")
#
#
# plt.figure()
# plt.imshow(viirs[0], vmin=0, vmax=60)
# plt.colorbar(label = "Intensite lumineuse [$nW/cm^2/sr$]")# plt.title('Viirs')
# plt.axis('off')
# plt.savefig('makeVIIRS/Results/viirs', dpi=300, bbox_inches="tight")
#
#
# plt.figure()
# plt.imshow(diff_viirs, vmin=0, vmax=60)
# plt.colorbar(label = "Intensite lumineuse [$nW/cm^2/sr$]")
# plt.axis('off')
# plt.savefig('makeVIIRS/Results/viirs_private', dpi=300, bbox_inches="tight")
#
#
# viirs[0][:] = diff_viirs
# viirs.save("diff_viirs")

plt.figure()
plt.imshow(sim_total, cmap='inferno')
plt.colorbar(label='Fraction public')
plt.axis('off')
plt.savefig('makeVIIRS/Results/viirs_public_total', bbox_inches='tight')
plt.close()
