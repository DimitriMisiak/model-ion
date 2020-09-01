#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:43:33 2020

@author: filippini
"""

import numpy as np
import matplotlib.pyplot as plt

e = 1.6e-19

def resolution(noise_f, Z, freq):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en frequentielle avec la méthode des trapèzes.

    Parameters
    ---------
    noise_f : np.array
        Bruit fréquentiel (PSD) du module d'amplification en :math:`V^2/Hz`

        Valeur de la fmax d'integration pour le calcul de la resolution

    Returns
    ----------
    res : float
        resolution du système d'amplification en :math:`eV`
    """
        
    noise_f = noise_f**2
    Z = np.abs(Z)
    signal_f = 333 * e * Z
    for i in [0,9,99]:
        print('frequency ', freq[i], ' FFT ', signal_f[i])
    NEPsquare2 = noise_f/(signal_f**2)
    # la case 0 correspond à 1Hz
    NEPsquare2 = 4/NEPsquare2
    res1 = np.sum(NEPsquare2) ** (-0.5)
    reso = res1*1e3  # On passe en eV
    return reso