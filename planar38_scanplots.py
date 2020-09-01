#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:38:58 2020

@author: misiak
"""


import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# =============================================================================
# Extra detector capacitance
# =============================================================================

capa_array = np.array([3, 8, 13, 28])
res_array = np.array([34.6, 32.9, 34, 42])

plt.figure()
plt.plot(capa_array, res_array, label='B', marker='o')

plt.ylabel('Resolution [eV]')
plt.xlabel('Hemt Line Equivalent Capacitance [pF]')
plt.grid()


# =============================================================================
# ALPHA multiplicative factor
# =============================================================================
alpha_array = np.array([0.25, 0.5, 0.75, 1, 1.5, 2])
res_array = np.array([19, 24, 29, 34, 44, 53])

plt.figure()
plt.plot(alpha_array, res_array, label='B', marker='o')

plt.ylabel('Resolution [eV]')
plt.xlabel('Capacitance matrix Multiplier Factor')

plt.ylabel('Resolution [eV]')
plt.xlabel('Multiplicative factor')
plt.grid()



# =============================================================================
# BETA multiplicative factor
# =============================================================================
alpha_array = np.array([0.25, 0.5, 0.75, 1, 1.5, 2])
res_array = np.array([20.6, 25, 29.5, 34, 43, 52.2])

plt.figure()
plt.plot(alpha_array, res_array, label='B', marker='o')

plt.ylabel('Resolution [eV]')
plt.xlabel('Capacitance matrix Multiplier Factor')

plt.ylabel('Resolution [eV]')
plt.xlabel('Multiplicative factor')
plt.grid()


