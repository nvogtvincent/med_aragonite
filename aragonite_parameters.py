#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 20:18:58 2021

@author: noam
"""
import numpy as np
import matplotlib.pyplot as plt

def fmt(x):
    s = f'{x:.1f}'
    if s.endswith('0'):
        s = f'{x:.0f}'
    return f'{s}'

T = np.array([5, 25, 37])

n = np.array([0.4, 1.7, 2.4])
nstdv = np.array([0.1, 0.1, 0.1])


# Note: (24/1000) conversion factor is to convert umol m-2 h-1 -> mmol m-2 d-1
k = np.array([21.8, 40.6, 45.1])*(24/1000)
kstdv = np.array([1.3, 1.1, 1.1])*(24/1000)

n_mean = np.polyfit(T, n, deg=1)
n_mean_p = np.poly1d(n_mean)
n_hi = np.polyfit(T, n+nstdv, deg=1)
n_hi_p = np.poly1d(n_hi)
n_lo = np.polyfit(T, n-nstdv, deg=1)
n_lo_p = np.poly1d(n_lo)

k_mean = np.polyfit(T, k, deg=2)
k_mean_p = np.poly1d(k_mean)
k_hi = np.polyfit(T, k+kstdv, deg=2)
k_hi_p = np.poly1d(k_hi)
k_lo = np.polyfit(T, k-kstdv, deg=2)
k_lo_p = np.poly1d(k_lo)

T_cont = np.linspace(0, 40, num=100)

# Plotting
plt.rcParams.update({'font.size' : 20})
f, ax = plt.subplots(2, 1, figsize=(15, 15))

ax[0].errorbar(T, n, yerr=nstdv, c='k', marker='.', linestyle='none', label='n (order)')
ax[0].plot(T_cont, n_mean_p(T_cont), 'k')
ax[0].plot(T_cont, n_hi_p(T_cont), 'k--')
ax[0].plot(T_cont, n_lo_p(T_cont), 'k--')
ax[0].set_xlabel('Temperature (C)')
ax[0].legend()

ax[1].errorbar(T, k, yerr=kstdv, c='b', marker='.', linestyle='none', label='k (rate constant, mmol/m2/d)')
ax[1].plot(T_cont, k_mean_p(T_cont), 'b')
ax[1].plot(T_cont, k_hi_p(T_cont), 'b--')
ax[1].plot(T_cont, k_lo_p(T_cont), 'b--')
ax[1].set_xlabel('Temperature (C)')
ax[1].legend()

plt.savefig('parameters.png', dpi=300)


# Now plot phase diagrams
Omega = np.linspace(1.5, 5, num=100)
T     = np.linspace(15, 35, num=100)
T, Omega = np.meshgrid(T, Omega)

R_lo = k_lo_p(T)*((Omega-1)**(n_lo_p(T)))
R_hi = k_hi_p(T)*((Omega-1)**(n_hi_p(T)))
R_mean = k_mean_p(T)*((Omega-1)**(n_mean_p(T)))

f, ax = plt.subplots(1, 1, figsize=(15, 15))
levels=np.linspace(0, 30, num=31, dtype=np.int8)

c_mean = ax.contour(T, Omega, R_mean, levels=levels, colors='k', linestyles='-', label='mean')
ax.clabel(c_mean, levels, fmt=fmt)

# c_lo = ax.contour(T, Omega, R_lo, levels=levels, colors='b', linestyles='-', label='low end')
# ax.clabel(c_lo, levels, fmt=fmt)

# c_hi = ax.contour(T, Omega, R_hi, levels=levels, colors='r', linestyles='-', label='high end')
# ax.clabel(c_hi, levels, fmt=fmt)

font = {'size': 20}
ax.set_xlabel('Temperature (C)')
ax.set_ylabel('Omega')

f.subplots_adjust(hspace=0.0, top=0.95)

f.suptitle('Inorganic aragonite precpitation rate (mmol m-2 d-1) ')
plt.savefig('/figures/parameter_space.png', dpi=300)


