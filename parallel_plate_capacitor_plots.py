#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:50:29 2020

@author: misiak
"""


import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

A = np.pi * (15e-3)**2

D = 10e-3

eps0 = 8.85e-12

epsr_air = 1
epsr_ge = 16.3

C = eps0 * epsr_ge * A / D

d_array = np.linspace(1e-3, 100e-3, int(1e4))
# a_array = np.linspace(np.pi * (1e-3)**2, np.pi * (100e-3)**2, int(1e4))
a_array = np.pi * (np.linspace(1e-3, 100e-3, int(1e4)))**2

capa_dist_ge = eps0 * epsr_ge * A / d_array
capa_dist_air = eps0 * epsr_air * A / d_array

capa_area_ge = eps0 * epsr_ge * a_array / D
capa_area_air = eps0 * epsr_air * a_array / D

fig, axes = plt.subplots(ncols=2, sharey=True)

ax = axes[0]
ax.plot(d_array, capa_dist_air, label='Air')
ax.plot(d_array, capa_dist_ge, label='Ge')
ax.plot(D, C, marker='o', ls='none', label='Ref')
ax.set_ylabel('Capacitance [F]')
ax.set_xlabel('Distance [m]')

ax = axes[1]
ax.plot(a_array, capa_area_air, label='Air')
ax.plot(a_array, capa_area_ge, label='Ge')
ax.plot(A, C, marker='o', ls='none', label='Ref')
ax.set_xlabel('Area [m${}^2$]')

for ax in axes:
    ax.grid()
    ax.legend()

fig.tight_layout()
fig.subplots_adjust(wspace=.0)