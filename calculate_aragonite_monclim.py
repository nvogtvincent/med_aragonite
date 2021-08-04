#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates a climatology from a decade of interest from the output
of calculate_aragonite.py
@author: Noam Vogt-Vincent
"""

import numpy as np
import os
from netCDF4 import Dataset
from datetime import datetime

##############################################################################
# Parameters #################################################################
##############################################################################

# File name (without extension)
fh = 'TA-pCO2_HISTORICAL'

# Decade of interest
dec = 1850       # First year of the decade (i.e. 1850 = 1850-1859)
clim_length = 50 # Length of climatology (years, default 10)

##############################################################################
# Parameter processing #######################################################
##############################################################################

root_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
out_fh = (root_dir + 'processed_data/' + fh +'_' +
          str(dec) + '-' + str(dec+clim_length-1) + '_MonClim.nc')
fh = root_dir + 'processed_data/' + fh + '.nc'

with Dataset(fh, mode='r') as nc:
    # Load variables
    lon = np.array(nc.variables['lon'][:])
    lat = np.array(nc.variables['lat'][:])
    time = np.array(nc.variables['time'][:]) # Days from Jan-1 1850

    # Extract time indices
    days_since_1850 = (dec - 1850)*360 + 15
    t0 = np.searchsorted(time, days_since_1850)
    t1 = t0 + 12*clim_length

    omega = nc.variables['omega'][:, t0:t1, :, :]
    rate = nc.variables['rate'][:, t0:t1, :, :]
    sst = nc.variables['sst'][:, t0:t1, :, :]

    fill = nc.variables['omega']._FillValue
    model_names = nc.model_names
    scenario    = nc.scenario
    par1        = nc.par1
    par2        = nc.par2

nmodel = np.shape(omega)[0]

# Calculate climatologies
omega_clim = np.zeros((nmodel, 12, len(lat), len(lon)))
rate_clim  = np.zeros_like(omega_clim)
sst_clim   = np.zeros_like(omega_clim)

coef = 1/(clim_length)

for t in range(12*clim_length):
    month = t%12
    omega_clim[:, month, :, :] += coef*omega[:, t, :, :]
    rate_clim[:, month, :, :] += coef*rate[:, t, :, :]
    sst_clim[:, month, :, :] += coef*sst[:, t, :, :]

# Mask variables
mask = np.ma.getmask(omega)[0, 0, :, :]
mask = np.tile(mask, (nmodel, 12, 1, 1))

omega_clim = np.ma.masked_array(omega_clim, mask=mask).filled(fill)
rate_clim = np.ma.masked_array(rate_clim, mask=mask).filled(fill)
sst_clim = np.ma.masked_array(sst_clim, mask=mask).filled(fill)

##############################################################################
# Write file #################################################################
##############################################################################

with Dataset(out_fh, mode='w') as nc:
    # Create dimensions
    nc.createDimension('lon', len(lon))
    nc.createDimension('lat', len(lat))
    nc.createDimension('month', 12)
    nc.createDimension('model', nmodel)

    # Create variables
    nc.createVariable('lon', 'f4', ('lon'), zlib=True)
    nc.createVariable('lat', 'f4', ('lat'), zlib=True)
    nc.createVariable('month', 'f4', ('month'), zlib=True)
    nc.createVariable('sst_clim', 'f4',
                      ('model', 'month', 'lat', 'lon'),
                      zlib=True, fill_value=fill)
    nc.createVariable('omega_clim', 'f4',
                      ('model', 'month', 'lat', 'lon'),
                      zlib=True, fill_value=fill)
    nc.createVariable('rate_clim', 'f4',
                      ('model', 'month', 'lat', 'lon'),
                      zlib=True, fill_value=fill)

    # Write non-model-dependent variables
    nc.variables['lon'].long_name = 'longitude'
    nc.variables['lon'].units = 'degrees_east'
    nc.variables['lon'].standard_name = 'longitude'
    nc.variables['lon'][:] = lon

    nc.variables['lat'].long_name = 'latitude'
    nc.variables['lat'].units = 'degrees_north'
    nc.variables['lat'].standard_name = 'latitude'
    nc.variables['lat'][:] = lat

    nc.variables['month'].long_name = 'month'
    nc.variables['month'].units = 'months'
    nc.variables['month'].standard_name = 'month'
    nc.variables['month'][:] = np.arange(1,13)

    nc.variables['sst_clim'].long_name = 'sea_surface_temperature_monthly_climatology'
    nc.variables['sst_clim'].units = 'degrees Celsius'
    nc.variables['sst_clim'].standard_name = 'sea_surface_temperature_climatology'
    nc.variables['sst_clim'][:] = sst_clim

    nc.variables['omega_clim'].long_name = 'aragonite_saturation_state_monthly_climatology'
    nc.variables['omega_clim'].units = 'no_units'
    nc.variables['omega_clim'].standard_name = 'omega_aragonite_climatology'
    nc.variables['omega_clim'][:] = omega_clim

    nc.variables['rate_clim'].long_name = 'aragonite_precipitation_rate_monthly_climatology'
    nc.variables['rate_clim'].units = 'mmol m-2 day-1'
    nc.variables['rate_clim'].standard_name = 'aragonite_precipitation_rate_climatology'
    nc.variables['rate_clim'][:] = rate_clim

    # Create global attributes
    date = datetime.now()
    date = date.strftime('%d/%m/%Y, %H:%M:%S')

    nc.model_names = model_names
    nc.scenario = scenario
    nc.date_created = date
    nc.par1 = par1
    nc.par2 = par2
    nc.decade = str(dec)


