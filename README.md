# lamps_conversions

Analyse lights conversions plans impact on light polution.

### Getting started
This repo is used to generate converions plans from the illumina model. For more inforrmation follow the Illumina user [guide](https://lx02.cegepsherbrooke.qc.ca/~aubema/index.php/Prof/Page).
In order to compare the impacts of differents conversions scenario, we suggest to first modelize the actual scenario.


### Obtaining private lights contributions
A challenge from modelizing light polutions comes from obtaining a private light inventory. We solve this problem by creating a light inventory from the VIIRS sattelites and we substract the public inventory.

Use the ```make_viirs``` scripts to simulate a fake viirs images from your lamps inventory.


### Simulating conversions scenario
You can create a new scenario by adding a file in the ```/screnarios``` folder.
The scripts conversions_plan.py calcul the sky intensity (zenith only) for the observer position.

### Ploting results
The following scripts are used to plot results from the modelisations
- ```plot_viirs.py```: Plot the simulate VIIRS image
- ```plot_all-sky.py``` : Generate a sky map from the observer positions
- ```plot_contrib-map.py``` : Generate a domain map with the pixel contribution to artificial radiance
- ```plot_spectral_bands``` : Plot the spectral bands of the lights simulate
