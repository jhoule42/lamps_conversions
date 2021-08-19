#!/usr/bin/env python3

import numpy as np
import MultiScaleData as MSD
from glob import glob
import yaml


with open('Inputs/lamps.lst') as f:
    lamps = f.read().split()


with open('domain.ini') as f:
    p = yaml.safe_load(f)
ymin = p['extents'][0]['ymin']
ymax = p['extents'][0]['ymax']
scale_min = np.deg2rad(p['scale_min']) # angle deg

arr = np.linspace(ymax, ymin, 84)
arr = (arr[1:] + arr[:-1]) / 2  # Trouver le point millieu
lat = np.zeros((83, 83))
lat[:] = arr[:, None]


R = 637813700  # Rayon de la terre (cm)
dx = scale_min * R * np.cos(np.deg2rad(lat))
dy = scale_min * R
S = dx*dy


with open("inputs_params.in") as f:
    p = yaml.safe_load(f)
alpha = p['angstrom_coefficient']
AOD = p['aerosol_optical_depth']
wav = np.loadtxt("Inputs/wav.lst")
refl = np.loadtxt('Inputs/refl.lst')

lims = np.loadtxt("Inputs/integration_limits.dat",skiprows=1)
wl,viirs = np.loadtxt('Lights/viirs.dat',skiprows=1).T
band = (wl>=lims[0]) & (wl<lims[1])
R = np.sum(viirs[band]) * (wl[1] - wl[0]) / (lims[1] - lims[0])

T_m = np.exp( -1 / ( (wav/1000)**4 * (115.6406 - 1.335/(wav/1000)**2) ) )
T_a = np.exp( -AOD * (wav/500)**(-alpha) )

altlp = MSD.Open("Inputs/PMB_altlp.hdf5")[0]
obstd = MSD.Open("Inputs/PMB_obstd.hdf5")[0]
obsth = MSD.Open("Inputs/PMB_obsth.hdf5")[0]
obstf = MSD.Open("Inputs/PMB_obstf.hdf5")[0]
phi_down = np.arctan2(obstd, obsth)
phi_up = np.arctan2(obstd, obsth-altlp)

viirs_lamps = []
for lamp in lamps:
    fctem,angle = np.loadtxt(f'Inputs/fctem_wl_690_lamp_{lamp}.dat').T
    sint = np.sin(np.deg2rad(angle))
    fctem /= 2*np.pi * np.sum( (fctem*sint) ) * np.deg2rad( angle[1] - angle[0] )

    T = (T_m*T_a)  ** (1 / np.cos(np.deg2rad(angle)))
    O_down = np.ones( angle.shape + obsth.shape ) - obstf
    O_down[ angle[:,None,None] < phi_down ] = 1
    O_down = O_down.transpose(1,2,0)
    O_up = np.ones( angle.shape + obsth.shape ) - obstf
    O_up[ angle[:,None,None] < phi_up ] = 1
    O_up = O_up.transpose(1,2,0)

    Gup = np.sum((O_up*fctem*T*sint)[:,:,angle<=70],2) / np.sum(sint[angle<=70])
    Fdown = 2*np.pi * np.sum( (fctem*sint)[angle>90] ) * np.deg2rad( angle[1] - angle[0])
    Fmoy = Fdown * np.sum((O_down*T*sint)[:,:,angle<=70],2) / np.sum(sint[angle<=70])

    A = 1/np.pi * refl * Fmoy + Gup

    lumlp = MSD.Open(f'Inputs/PMB_690_lumlp_{lamp}.hdf5')
    phi = lumlp[0]*(lims[1]-lims[0]) * 1e9
    DNB = phi/S * R * A

    # Somme sur les lampes
    viirs_lamps.append(DNB)


sim_viirs = np.sum(viirs_lamps,0)
np.save('sim_viirs', sim_viirs)
