#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 14:37:13 2021
This script produces a netcdf file with aragonite saturation state, sst, and
aragonite precipitation rate from gridded input. A netcdf file is produced
for each scenario provided.
@author: Noam Vogt-Vincent
"""

import PyCO2SYS as pyco2
import numpy as np
import os
import seawater as sw
from netCDF4 import Dataset
from progress.bar import Bar
from datetime import datetime

##############################################################################
# Parameters #################################################################
##############################################################################

# Array of model names
model_names = ['UKESM1-0-LL_r1i1p1f2',
               'GFDL-CM4_r1i1p1f1',
               'CanESM5_r1i1p1f1']

# Array of scenario names
scenario_names = ['HISTORICAL',
                  'SSP245',
                  'SSP585']

# Variables
sst_var   = 'tos'       # Name of sea surface temperature variable
sss_var   = 'sos'       # Name of sea surface salinity variable

si_var    = 'sios'      # Name of the silicate variable
po4_var   = 'po4'       # Name of the phosphate variable
use_nutri = False        # Whether to use Si/PO4 in calculations

par1_var  = 'talk'      # Name of the first carbon system parameter
par1_type = 1           # Parameter type
par2_var  = 'spco2'    # Name of the first carbon system parameter
par2_type = 4           # Parameter type

output_name = 'TA-pCO2'  # Output name

# Parameters for the rate parameterisation
Ki = lambda t : ((-0.0177*(t**2)) + (1.4697*t) + 14.893)*(24/1000)
Ni = lambda t : (0.0628*t + 0.0985)

##############################################################################
# Parameter processing #######################################################
##############################################################################

root_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
data_dir = root_dir + 'data/'
proc_dir = root_dir + 'processed_data/'

# Parameters for the rate parameterisation (derive from results in Burton &
# Walter, 1987)
T = np.array([5, 25, 37])
k_obs = np.array([21.8, 40.6, 45.1])*(24/1000) # umol/m2/h -> mmol/m2/d
n_obs = np.array([0.4, 1.7, 2.4])
k_int = np.polyfit(T, k_obs, deg=2)
n_int = np.polyfit(T, n_obs, deg=1)

Ki = lambda t : ((k_int[0]*(t**2)) + (k_int[1]*t) + k_int[2])
Ni = lambda t : (n_int[0]*t + n_int[1])

fh = lambda v, m, s : data_dir + m + '/' + v + '_' + m + '_' + s + '.nc'
fh_proc = lambda v, s : proc_dir + v + '_' + s + '.nc'

nmodel = len(model_names)
nscenario = len(scenario_names)

# Check that variable types are allowed!
if not (par1_type in [1, 2, 3, 4]):
    raise Exception('Parameter 1 type not allowed!')
elif not (par2_type in [1, 2, 3, 4]):
    raise Exception('Parameter 2 type not allowed!')

##############################################################################
# Load data ##################################################################
##############################################################################

print('Loading data...')

for i, scenario in enumerate(scenario_names):
    print('Scenario ' + str(i+1) + '/' + str(nscenario))
    out_fh = fh_proc(output_name, scenario)

    for j, model in enumerate(model_names):
        print('Model ' + str(j+1) + '/' + str(nmodel))

        sst_fh  = fh(sst_var, model, scenario)
        sss_fh  = fh(sss_var, model, scenario)
        par1_fh = fh(par1_var, model, scenario)
        par2_fh = fh(par2_var, model, scenario)

        with Dataset(sst_fh, mode='r') as nc:
            # If this is the first model, obtain lat/lon/time/mask
            if j == 0:
                lon = np.array(nc.variables['lon'][:])
                lat = np.array(nc.variables['lat'][:])
                time = np.array(nc.variables['time'][:]) # Days from Jan-1 1850
                time_bnds = np.array(nc.variables['time_bnds'][:])

                ntime = len(time) # Record number of time steps

                fill = nc.variables[sst_var]._FillValue
                mask = np.ma.getmask(nc.variables[sst_var][0, :, :])

            sst = np.squeeze(nc.variables[sst_var][:])

        with Dataset(sss_fh, mode='r') as nc:
            sss = np.squeeze(nc.variables[sss_var][:])

        with Dataset(par1_fh, mode='r') as nc:
            par1 = np.squeeze(nc.variables[par1_var][:])

        with Dataset(par2_fh, mode='r') as nc:
            par2 = np.squeeze(nc.variables[par2_var][:])

        if use_nutri:
            si_fh  = fh(si_var, model, scenario)
            po4_fh = fh(po4_var, model, scenario)

            with Dataset(si_fh, mode='r') as nc:
                si = np.squeeze(nc.variables[si_var][:])

            with Dataset(po4_fh, mode='r') as nc:
                po4 = np.squeeze(nc.variables[po4_var][:])

        # Further processing may be needed depending on the parameter type:
        # par_type = 1 (TA):   mol m-3 -> umol kg-1
        # par_type = 2 (DIC):  mol m-3 -> umol kg-1
        # par_type = 3 (pH):   no processing needed
        # par_type = 4 (pCO2): Pa -> uatm

        if ((1 in [par1_type, par2_type]) or (2 in [par1_type, par2_type])):
            # Calculate salinity if required for conversion
            rho = sw.dens0(sss, sst)
        elif use_nutri:
            rho = sw.dens0(sss, sst)

        if ((par1_type == 1) or (par1_type == 2)):
            # Convert mol m-3 -> umol kg-1
            par1 *= (1E6/rho)
        elif par1_type == 4:
            # Convert Pa -> uatm
            par1 *= (1E6/101325)

        if ((par2_type == 1) or (par2_type == 2)):
            # Convert mol m-3 -> umol kg-1
            par2 *= (1E6/rho)
        elif par2_type == 4:
            # Convert Pa -> uatm
            par2 *= (1E6/101325)

        if use_nutri:
            # Convert mol m-3 -> umol kg-1
            si  *= (1E6/rho)
            po4 *= (1E6/rho)

        # Now generate an empty array for omega
        omega = np.zeros_like(sst)

        # Calculate omega
        bar = Bar('Calculating aragonite saturation state...', max=len(time))
        for t in range(len(time)):
            if use_nutri:
                omega[t, :, :] = pyco2.sys(par1=par1[t, :, :],
                                           par1_type=par1_type,
                                           par2=par2[t, :, :],
                                           par2_type=par2_type,
                                           salinity=sss[t, :, :],
                                           temperature=sst[t, :, :],
                                           total_silicate=si[t, :, :],
                                           total_phosphate=po4[t, :, :]
                                           )['saturation_aragonite']
            else:
                omega[t, :, :] = pyco2.sys(par1=par1[t, :, :],
                                           par1_type=par1_type,
                                           par2=par2[t, :, :],
                                           par2_type=par2_type,
                                           salinity=sss[t, :, :],
                                           temperature=sst[t, :, :],
                                           )['saturation_aragonite']
            bar.next()
        bar.finish()

        # Generate the precipitation rate of aragonite
        rate = Ki(sst)*((omega-1)**Ni(sst))

        # Apply mask properly
        omega[np.isnan(omega)] = fill
        rate[np.isnan(rate)]   = fill
        rate[omega <= 1]       = fill # Since this is undefined
        sst                    = sst.filled(fill)

        # Save to netcdf
        print('Writing to file...')
        if j == 0:
            # If this is the first model, create the required dims/vars
            with Dataset(out_fh, mode='w') as nc:
                # Create dimensions
                nc.createDimension('lon', len(lon))
                nc.createDimension('lat', len(lat))
                nc.createDimension('time', len(time))
                nc.createDimension('bnds', 2)
                nc.createDimension('model', nmodel)

                # Create variables
                nc.createVariable('lon', 'f4', ('lon'), zlib=True)
                nc.createVariable('lat', 'f4', ('lat'), zlib=True)
                nc.createVariable('time', 'f4', ('time'), zlib=True)
                nc.createVariable('time_bnds', 'f4', ('time', 'bnds'),
                                  zlib=True)
                nc.createVariable('sst', 'f4',
                                  ('model', 'time', 'lat', 'lon'),
                                  zlib=True, fill_value=fill)
                nc.createVariable('omega', 'f4',
                                  ('model', 'time', 'lat', 'lon'),
                                  zlib=True, fill_value=fill)
                nc.createVariable('rate', 'f4',
                                  ('model', 'time', 'lat', 'lon'),
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

                nc.variables['time'].long_name = 'time'
                nc.variables['time'].units = 'days since 1850-01-01'
                nc.variables['time'].standard_name = 'time'
                nc.variables['time'].calendar = "360_day"
                nc.variables['time'].bounds = 'time_bnds'
                nc.variables['time'][:] = time

                nc.variables['time_bnds'].long_name = 'time_bounds'
                nc.variables['time_bnds'].units = 'days since 1850-01-01'
                nc.variables['time_bnds'].standard_name = 'time_bounds'
                nc.variables['time_bnds'][:] = time_bnds

                nc.variables['sst'].long_name = 'sea_surface_temperature'
                nc.variables['sst'].units = 'degrees Celsius'
                nc.variables['sst'].standard_name = 'sea_surface_temperature'

                nc.variables['omega'].long_name = 'aragonite_saturation_state'
                nc.variables['omega'].units = 'no_units'
                nc.variables['omega'].standard_name = 'omega_aragonite'

                nc.variables['rate'].long_name = 'aragonite_precipitation_rate'
                nc.variables['rate'].units = 'mmol m-2 day-1'
                nc.variables['rate'].standard_name = 'aragonite_precipitation_rate'

                # Create global attributes
                namelist = ''
                for k, model in enumerate(model_names):
                    namelist = namelist + str(k+1) + ':' + model + '/'

                date = datetime.now()
                date = date.strftime('%d/%m/%Y, %H:%M:%S')

                nc.model_names = namelist
                nc.scenario = scenario
                nc.date_created = date
                nc.par1 = par1_var
                nc.par2 = par2_var

        with Dataset(out_fh, mode='r+') as nc:
            # Write model variables
            nc.variables['sst'][j, :, :, :] = sst[:ntime, :, :]
            nc.variables['omega'][j, :, :, :] = omega[:ntime, :, :]
            nc.variables['rate'][j, :, :, :] = rate[:ntime, :, :]

print('Complete!')
