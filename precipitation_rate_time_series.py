#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 16:29:16 2021
This script calculates the mean, minimum and maximum precipitation rate for
a region, and separates the omega and sst contributions
@author: Noam Vogt-Vincent
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter1d

##############################################################################
# File locations #############################################################
##############################################################################

model_names = ['UKESM1-0-LL',
               'GFDL-CM4',
               'CanESM5']

root_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = root_dir + '/processed_data/'

fh = ['TA-pCO2_HISTORICAL',
      'TA-pCO2_SSP245',
      'TA-pCO2_SSP585']

which_scenario = 'SSP245'
fig_fh = root_dir + '/figures/Eastern_Med_PR_' + which_scenario + '.png'

fh_clim = fh.copy()
for i in range(len(fh)):
    fh_clim[i] = data_dir + fh[i] + '_1850-1899_MonClim.nc'
    fh[i] = data_dir + fh[i] + '.nc'



##############################################################################
# Parameters #################################################################
##############################################################################

lon_bnds = [32, 40]
lat_bnds = [30, 40]

nmodel = 3               # Number of models in file

# Parameters for the rate parameterisation (derive from results in Burton &
# Walter, 1987)
T = np.array([5, 25, 37])
k_obs = np.array([21.8, 40.6, 45.1])*(24/1000) # umol/m2/h -> mmol/m2/d
n_obs = np.array([0.4, 1.7, 2.4])
k_int = np.polyfit(T, k_obs, deg=2)
n_int = np.polyfit(T, n_obs, deg=1)


Ki = lambda t : ((k_int[0]*(t**2)) + (k_int[1]*t) + k_int[2])
Ni = lambda t : (n_int[0]*t + n_int[1])

##############################################################################
# Load data ##################################################################
##############################################################################

print('Setting up...')
with Dataset(fh[0], mode='r') as nc:
    lon = np.array(nc.variables['lon'][:])
    lat = np.array(nc.variables['lat'][:])
    time_hist = np.array(nc.variables['time'][:]) # Days since Jan-1 1850

    ntime_hist = len(time_hist)
    nlon = len(lon)
    nlat = len(lat)

with Dataset(fh[1], mode='r') as nc:
    time_future = np.array(nc.variables['time'][:]) # Days since Jan-1 1850

    ntime_future = len(time_future)

time = np.append(time_hist, time_future)
ntime = ntime_hist + ntime_future

# Import the monthly climatology
with Dataset(fh_clim[0], mode='r') as nc:
    sst_clim = nc.variables['sst_clim'][:]
    rate_clim = nc.variables['rate_clim'][:]
    omega_clim = nc.variables['omega_clim'][:]

##############################################################################
# Process data ###############################################################
##############################################################################

# Calculate a mask for the region of interest
lon, lat = np.meshgrid(lon, lat)

mask0 = ~((lon > lon_bnds[0])*
                (lon < lon_bnds[1])*
                (lat > lat_bnds[0])*
                (lat < lat_bnds[1]))

# Calculate a monthly climatology
sst_clim = np.zeros((3, 12, nlat, nlon))
omega_clim = np.zeros_like(sst_clim)

# Create a new array with the following information:
# int_rate(model, assumption, time)
# where assumptions = (0) var omega/sst (rate), (1) fixed sst, (2) fixed omega
pr = np.zeros((nmodel, 3, ntime))

# HISTORICAL
for model in range(nmodel):
    for assumption in range(3):
        if assumption == 0:
            with Dataset(fh[0], mode='r') as nc:
                rate_ = nc.variables['rate'][model, :, :, :]

            rate_ = np.ma.masked_array(rate_, mask=np.tile(mask0,
                                                           (ntime_hist, 1, 1)))
            pr[model, assumption, :ntime_hist] = np.mean(rate_, (1,2))
        elif assumption == 1:
            # Varying omega, climatological SST
            with Dataset(fh[0], mode='r') as nc:
                omega_ = nc.variables['omega'][model, :, :, :]
                omega_ = np.ma.masked_array(omega_, mask=np.tile(mask0,
                                                           (ntime_hist, 1, 1)))
            with Dataset(fh_clim[0], mode='r') as nc:
                sst_ = nc.variables['sst_clim'][model, :, :, :]
                sst_ = np.tile(sst_, (int(ntime_hist/12), 1, 1))
                sst_ = np.ma.masked_array(sst_, mask=np.tile(mask0,
                                                           (ntime_hist, 1, 1)))

            rate_ = Ki(sst_)*(omega_-1)**(Ni(sst_))
            pr[model, assumption, :ntime_hist] = np.mean(rate_, (1,2))
        elif assumption == 2:
            # Varying SST, climatological omega
            with Dataset(fh[0], mode='r') as nc:
                sst_ = nc.variables['sst'][model, :, :, :]
                sst_ = np.ma.masked_array(sst_, mask=np.tile(mask0,
                                                             (ntime_hist, 1, 1)))
            with Dataset(fh_clim[0], mode='r') as nc:
                omega_ = nc.variables['omega_clim'][model, :, :, :]
                omega_ = np.tile(omega_, (int(ntime_hist/12), 1, 1))
                omega_ = np.ma.masked_array(omega_, mask=np.tile(mask0,
                                                           (ntime_hist, 1, 1)))

            rate_ = Ki(sst_)*(omega_-1)**(Ni(sst_))
            pr[model, assumption, :ntime_hist] = np.mean(rate_, (1,2))

# FUTURE
if which_scenario == 'SSP245':
    which_scenario_num = 1
elif which_scenario == 'SSP585':
    which_scenario_num = 2
else:
    print('Oh dear, scenario not understood!')

for model in range(nmodel):
    for assumption in range(3):
        if assumption == 0:
            with Dataset(fh[which_scenario_num], mode='r') as nc:
                rate_ = nc.variables['rate'][model, :, :, :]

            rate_ = np.ma.masked_array(rate_, mask=np.tile(mask0,
                                                           (ntime_future, 1, 1)))
            pr[model, assumption, ntime_hist:] = np.mean(rate_, (1,2))
        elif assumption == 1:
            # Varying omega, climatological SST
            with Dataset(fh[which_scenario_num], mode='r') as nc:
                omega_ = nc.variables['omega'][model, :, :, :]
                omega_ = np.ma.masked_array(omega_, mask=np.tile(mask0,
                                                           (ntime_future, 1, 1)))
            with Dataset(fh_clim[0], mode='r') as nc:
                sst_ = nc.variables['sst_clim'][model, :, :, :]
                sst_ = np.tile(sst_, (int(ntime_future/12), 1, 1))
                sst_ = np.ma.masked_array(sst_, mask=np.tile(mask0,
                                                           (ntime_future, 1, 1)))

            rate_ = Ki(sst_)*(omega_-1)**(Ni(sst_))
            pr[model, assumption, ntime_hist:] = np.mean(rate_, (1,2))
        elif assumption == 2:
            # Varying SST, climatological omega
            with Dataset(fh[which_scenario_num], mode='r') as nc:
                sst_ = nc.variables['sst'][model, :, :, :]
                sst_ = np.ma.masked_array(sst_, mask=np.tile(mask0,
                                                             (ntime_future, 1, 1)))
            with Dataset(fh_clim[0], mode='r') as nc:
                omega_ = nc.variables['omega_clim'][model, :, :, :]
                omega_ = np.tile(omega_, (int(ntime_future/12), 1, 1))
                omega_ = np.ma.masked_array(omega_, mask=np.tile(mask0,
                                                           (ntime_future, 1, 1)))

            rate_ = Ki(sst_)*(omega_-1)**(Ni(sst_))
            pr[model, assumption, ntime_hist:] = np.mean(rate_, (1,2))

# Create a new array with the following information:
# int_rate(model, assumption, range, time)
# where assumptions = (0) var omega/sst (rate), (1) fixed sst, (2) fixed omega
# where range = (0) min, (1) mean, (2) max
tyear = 1850 + np.linspace(0, 250, num=251)
nyear = len(tyear)

pr_ann = np.zeros((nmodel, 3, 3, nyear))
for model in range(nmodel):
    for t in range(nyear):
        for i in range(3):
            pr_ann[model, :, 0, t] = np.amin(pr[model, :, (t*12):((t+1)*12)],
                                             axis=1)
            pr_ann[model, :, 1, t] = np.mean(pr[model, :, (t*12):((t+1)*12)],
                                              axis=1)
            pr_ann[model, :, 2, t] = np.amax(pr[model, :, (t*12):((t+1)*12)],
                                             axis=1)

pr_ann = gaussian_filter1d(pr_ann, 5)
##############################################################################
# Plotting ###################################################################
##############################################################################
plt.rcParams.update({'font.size' : 10})

tyear = 1850 + np.linspace(0, 250, num=251)
f, ax = plt.subplots(3, 1, figsize=(15, 15))

for i in range(nmodel):

    ax[i].fill_between(tyear,
                    pr_ann[i, 2, 0, :],
                    pr_ann[i, 2, 2, :],
                    facecolor='r',
                    alpha=0.2)
    ax[i].plot(tyear, pr_ann[i, 2, 1, :], 'r-', linewidth=0.5,
            label='Varying SST, PI climatological Omega')

    ax[i].fill_between(tyear,
                    pr_ann[i, 1, 0, :],
                    pr_ann[i, 1, 2, :],
                    facecolor='b',
                    alpha=0.2)
    ax[i].plot(tyear, pr_ann[i, 1, 1, :], 'b-', linewidth=0.5,
            label='Varying Omega, PI climatological SST')

    ax[i].fill_between(tyear,
                    pr_ann[i, 0, 0, :],
                    pr_ann[i, 0, 2, :],
                    facecolor='k',
                    alpha=0.2)
    ax[i].plot(tyear, pr_ann[i, 0, 1, :], 'k-', linewidth=0.5,
            label='Varying Omega and SST')

    if i == 2:
        ax[i].set_xlabel('Year')
    else:
        ax[i].set_xticks([])

    ax[i].set_ylabel('Precipitation rate (mmol/m2/d)')

    ax[i].set_xlim(1900, 2100)
    ax[i].title.set_text(model_names[i])

title = f.suptitle('Aragonite precipitation rate under CMIP6 Historical and ScenarioMIP ' + which_scenario + ' in the eastern Mediterranean')
f.subplots_adjust(hspace=0.0, top=0.95)
# plt.legend()

plt.savefig(fig_fh, dpi=300)
