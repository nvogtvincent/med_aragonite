#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 19:14:41 2021
This script plots CMIP6 TA and pH against measurements
@author: Noam Vogt-Vincentg
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset

##############################################################################
# File locations #############################################################
##############################################################################
root_dir = os.path.dirname(os.path.realpath(__file__)) + '/..'
data_dir = root_dir + '/data/'

model_names = ['UKESM1-0-LL_r1i1p1f2',
               'GFDL-CM4_r1i1p1f1',
               'CanESM5_r1i1p1f1']

fh_talk = [data_dir + model_names[0] + '/talk_UKESM1-0-LL_r1i1p1f2_HISTORICAL.nc',
           data_dir + model_names[1] + '/talk_GFDL-CM4_r1i1p1f1_HISTORICAL.nc',
           data_dir + model_names[2] + '/talk_CanESM5_r1i1p1f1_HISTORICAL.nc']
fh_phos = [data_dir + model_names[0] + '/phos_UKESM1-0-LL_r1i1p1f2_HISTORICAL.nc',
           data_dir + model_names[1] + '/ph_GFDL-CM4_r1i1p1f1_HISTORICAL.nc',
           data_dir + model_names[1] + '/ph_GFDL-CM4_r1i1p1f1_HISTORICAL.nc']   # Not used/plotted, using GFDL as placeholder to avoid crash
fh_sst  = [data_dir + model_names[0] + '/tos_UKESM1-0-LL_r1i1p1f2_HISTORICAL.nc',
           data_dir + model_names[1] + '/tos_GFDL-CM4_r1i1p1f1_HISTORICAL.nc',
           data_dir + model_names[2] + '/tos_CanESM5_r1i1p1f1_HISTORICAL.nc']

fig_fh = root_dir + '/figures/CMIP6_obs_TA_pH.png'

obs_fh = root_dir + '/data/obs/EMed_pH_TA.csv'
era_fh = root_dir + '/data/ERA5/ERA5_SST.nc'

nmodel = len(model_names)

##############################################################################
# Parameters #################################################################
##############################################################################

lon_bnds = [32, 40]
lat_bnds = [30, 40]

# Density for unit conversion
rho0 = 1027.2

##############################################################################
# Load data ##################################################################
##############################################################################

# Load the observations
raw_obs = np.genfromtxt(obs_fh, delimiter=',', skip_header=True)

t_obs = raw_obs[:, 2]
pH_obs = raw_obs[:, 4]
TA_obs = raw_obs[:, 5]

t_obs_d13C = t_obs[:119]
t_obs_dir1 = t_obs[119:152]
t_obs_dir2 = t_obs[152:]

pH_obs_d13C = pH_obs[:119]
pH_obs_dir1 = pH_obs[119:152]
pH_obs_dir2 = pH_obs[152:]

TA_obs_d13C = TA_obs[:119]
TA_obs_dir1 = TA_obs[119:152]
TA_obs_dir2 = TA_obs[152:]

with Dataset(era_fh, mode='r') as nc:
    lon_era = nc.variables['longitude'][:]
    lat_era = nc.variables['latitude'][:]
    sst_era = nc.variables['sst'][:504,0,:,:] - 273.13
    time_era = nc.variables['time'][:504]
    ntime_era = len(time_era)

lon_era, lat_era = np.meshgrid(lon_era, lat_era)
mask_era = ~((lon_era > lon_bnds[0])*
                (lon_era < lon_bnds[1])*
                (lat_era > lat_bnds[0])*
                (lat_era < lat_bnds[1]))

sst_era = np.ma.masked_array(sst_era, mask=np.tile(mask_era,
                                                   (ntime_era, 1, 1)))
sst_era = np.mean(sst_era, (1,2))
tyear_era = np.linspace(1979 + 0.5, 2021 - 0.5, num=42)
nyear_era = len(tyear_era)

sst_era_ann = np.zeros((3, nyear_era))

for t in range(nyear_era):
    sst_era_ann[0, t] = np.min(sst_era[(t*12):((t+1)*12)])
    sst_era_ann[1, t] = np.mean(sst_era[(t*12):((t+1)*12)])
    sst_era_ann[2, t] = np.max(sst_era[(t*12):((t+1)*12)])

# Load the model
with Dataset(fh_talk[0], mode='r') as nc:
    lon = np.array(nc.variables['lon'][:])
    lat = np.array(nc.variables['lat'][:])
    time = np.array(nc.variables['time'][:]) # Days since Jan-1 1850

    ntime = len(time)
    nlon = len(lon)
    nlat = len(lat)

# Calculate a mask for the region of interest
lon, lat = np.meshgrid(lon, lat)

mask0 = ~((lon > lon_bnds[0])*
                (lon < lon_bnds[1])*
                (lat > lat_bnds[0])*
                (lat < lat_bnds[1]))

# Create new arrays with the following information:
# pH(model, time)
# TA(model, time)
# SST(model, time)
pH_model = np.zeros((nmodel, ntime))
TA_model = np.zeros((nmodel, ntime))
SST_model = np.zeros((nmodel, ntime))

for model in range(nmodel):
    with Dataset(fh_talk[model], mode='r') as nc:
        TA = np.squeeze(nc.variables['talk'][:])
        TA *= (1E6/rho0)

    try:
        with Dataset(fh_phos[model], mode='r') as nc:
            pH = np.squeeze(nc.variables['phos'][:])
    except:
        with Dataset(fh_phos[model], mode='r') as nc:
            pH = np.squeeze(nc.variables['ph'][:])

    with Dataset(fh_sst[model], mode='r') as nc:
        SST = np.squeeze(nc.variables['tos'][:])


    TA = np.ma.masked_array(TA, mask=np.tile(mask0,
                                             (ntime, 1, 1)))
    pH = np.ma.masked_array(pH, mask=np.tile(mask0,
                                             (ntime, 1, 1)))
    SST = np.ma.masked_array(SST, mask=np.tile(mask0,
                                               (ntime, 1, 1)))
    TA_model[model, :] = np.mean(TA, (1,2))
    pH_model[model, :] = np.mean(pH, (1,2))
    SST_model[model, :] = np.mean(SST, (1,2))

# Now calculate annual range
tyear = np.linspace(1850 + 0.5, 2015 - 0.5, num=165)
nyear = len(tyear)

pH_model_ann = np.zeros((nmodel, 3, nyear))
TA_model_ann = np.zeros((nmodel, 3, nyear))
SST_model_ann = np.zeros((nmodel, 3, nyear))

for model in range(nmodel):
    for t in range(nyear):
        pH_model_ann[model, 0, t] = np.min(pH_model[model, (t*12):((t+1)*12)])
        pH_model_ann[model, 1, t] = np.mean(pH_model[model, (t*12):((t+1)*12)])
        pH_model_ann[model, 2, t] = np.max(pH_model[model, (t*12):((t+1)*12)])

        TA_model_ann[model, 0, t] = np.min(TA_model[model, (t*12):((t+1)*12)])
        TA_model_ann[model, 1, t] = np.mean(TA_model[model, (t*12):((t+1)*12)])
        TA_model_ann[model, 2, t] = np.max(TA_model[model, (t*12):((t+1)*12)])

        SST_model_ann[model, 0, t] = np.min(SST_model[model, (t*12):((t+1)*12)])
        SST_model_ann[model, 1, t] = np.mean(SST_model[model, (t*12):((t+1)*12)])
        SST_model_ann[model, 2, t] = np.max(SST_model[model, (t*12):((t+1)*12)])

# Generate a time axis for the model data
t_model = np.linspace(1850 + (1/24), 2015 - (1/24), num=1980)


##############################################################################
# Plotting ###################################################################
##############################################################################

f, ax = plt.subplots(3, 1, figsize=(15, 15))
f.suptitle('Proxy vs CMIP6 TA and pH comparison')
f.subplots_adjust(hspace=0.0, top=0.95)

ax[0].scatter(t_obs_d13C, pH_obs_d13C, c='k', s=20, marker='s', label='Proxy reconstructed pH')
ax[0].scatter(t_obs_dir1, pH_obs_dir1, c='g', s=20, marker='v', label='Directly measured pH')
ax[0].scatter(t_obs_dir2, pH_obs_dir2, c='r', s=20, marker='^', label='DIC derived pH')
ax[0].set_ylabel('Total pH')
ax[0].fill_between(tyear,
                    pH_model_ann[0, 0, :],
                    pH_model_ann[0, 2, :],
                    facecolor='r',
                    alpha=0.1)
ax[0].plot(tyear, pH_model_ann[0, 1, :], 'r-', linewidth=0.5,
            label='UKESM1-0-LL eastern Mediterranean mean surface pH')
ax[0].fill_between(tyear,
                    pH_model_ann[1, 0, :],
                    pH_model_ann[1, 2, :],
                    facecolor='b',
                    alpha=0.1)
ax[0].plot(tyear, pH_model_ann[1, 1, :], 'b-', linewidth=0.5,
            label='GFDL-CM4 eastern Mediterranean mean surface pH')
ax[0].legend()
ax[0].set_xlim(1945, 2016)

ax[1].scatter(t_obs_d13C, TA_obs_d13C, c='k', s=20, marker='s',  label='Proxy reconstructed TA')
ax[1].scatter(t_obs_dir1, TA_obs_dir1, c='g', s=20, marker='v',  label='pH/DIC-derived TA')
ax[1].scatter(t_obs_dir2, TA_obs_dir2, c='r', s=20, marker='^',  label='Directly measured TA')
ax[1].set_ylabel('Total Alkalinity (umol Kg-1)')
ax[1].set_xlabel('Year')
ax[1].fill_between(tyear,
                    TA_model_ann[0, 0, :],
                    TA_model_ann[0, 2, :],
                    facecolor='r',
                    alpha=0.1)
ax[1].plot(tyear, TA_model_ann[0, 1, :], 'r-', linewidth=0.5,
            label='UKESM1-0-LL eastern Mediterranean mean surface TA')
ax[1].fill_between(tyear,
                    TA_model_ann[1, 0, :],
                    TA_model_ann[1, 2, :],
                    facecolor='b',
                    alpha=0.1)
ax[1].plot(tyear, TA_model_ann[1, 1, :], 'b-', linewidth=0.5,
            label='GFDL-CM4 eastern Mediterranean mean surface TA')
ax[1].fill_between(tyear,
                    TA_model_ann[2, 0, :],
                    TA_model_ann[2, 2, :],
                    facecolor='g',
                    alpha=0.1)
ax[1].plot(tyear, TA_model_ann[2, 1, :], 'g-', linewidth=0.5,
            label='CanESM5 eastern Mediterranean mean surface TA')
ax[1].legend()
ax[1].set_xlim(1945, 2016)


ax[2].set_ylabel('SST (C)')
ax[2].set_xlabel('Year')
ax[2].fill_between(tyear,
                    SST_model_ann[0, 0, :],
                    SST_model_ann[0, 2, :],
                    facecolor='r',
                    alpha=0.1)
ax[2].plot(tyear, SST_model_ann[0, 1, :], 'r-', linewidth=0.5,
            label='UKESM1-0-LL eastern Mediterranean mean surface SST')
ax[2].fill_between(tyear,
                    SST_model_ann[1, 0, :],
                    SST_model_ann[1, 2, :],
                    facecolor='b',
                    alpha=0.1)
ax[2].plot(tyear, SST_model_ann[1, 1, :], 'b-', linewidth=0.5,
            label='GFDL-CM4 eastern Mediterranean mean surface SST')
ax[2].fill_between(tyear,
                    SST_model_ann[2, 0, :],
                    SST_model_ann[2, 2, :],
                    facecolor='g',
                    alpha=0.1)
ax[2].plot(tyear, SST_model_ann[2, 1, :], 'g-', linewidth=0.5,
            label='CanESM5 eastern Mediterranean mean surface SST')
ax[2].fill_between(tyear_era,
                    sst_era_ann[0, :],
                    sst_era_ann[2, :],
                    facecolor='k',
                    alpha=0.1)
ax[2].plot(tyear_era, sst_era_ann[1, :], 'k-', linewidth=0.5,
            label='ERA5 eastern Mediterranean mean surface SST')
ax[2].legend()
ax[2].set_xlim(1945, 2016)

plt.savefig(fig_fh, dpi=300)

