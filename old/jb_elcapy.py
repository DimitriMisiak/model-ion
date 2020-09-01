#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 13:01:29 2020

@author: misiak

Modelization of a planar detector with 2 electrodes.
Simplistic modelization taking into account only the capacitance from each
electronic line.
The objective is to evaluate how the noise coming from the second electrode
affects the first electrode measure.
"""

import lcapy
import numpy as np
import matplotlib.pyplot as plt

# plt.close('all')

# a = Circuit("""
#     R 1 0; down
#     W 1 2; right
#     C 2 0_2; down
#     W 0 0_2; right
# """)

# a.draw(style='european')

# b = a.noisy()

# b.draw(style='european')

# Vn = b.C.V.n

# Vns = Vn.subs({'R':10e3, 'C':100e-9, 'T':293, 'k':1.38e-23})



# vf = np.linspace(0, 10e3, 200)

# Vns_array = Vns(f).evaluate(vf)

# A = Vns(f)

# A.plot(vf, ylabel='ASD (nV/rootHz)')

# plt.figure()
# plt.plot(vf, Vns_array)
# plt.ylabel('ASD (nV/rootHz)')
# plt.xlabel('Fraquency [Hz]')
# plt.xscale('log')
# plt.yscale('log')
# plt.grid()


def Z_c(freq, C):
    return 1/(2*np.pi*freq*C*1j)

Cd = 10e-12
Cp = 10e-12
Rpolar = 1e10
Cc = 2e-9
Chemt = 5e-12


frequency = np.arange(1,100000,1)
Zd = Z_c(frequency, Cd)
Zp = Z_c(frequency, Cp)
Zc = Z_c(frequency, Cc)
Zhemt = Z_c(frequency, Chemt)

circuit = lcapy.Circuit("""
                        Cd 1 0; down
                        W 1 1_2; right
                        Cp 1_2 0_1 10e-12; down
                        W 1_2 1_3; right
                        Rpolar 1_3 0_3 1e10; down
                        Cc 1_3 2 2e-9; right
                        W 2 2_1; right
                        Chemt 2_1 0_2 5e-12; down
                        W 2_1 2_2; right
                        Vn 2_2 3 noise; right
                        In 2_2 0_4 noise; down
                        W 0 0_1; right
                        W 0_1 0_2; right
                        W 0_2 0_3; right
                        W 0_3 0_4; right
                        """)

# circuit = (((lcapy.C(Cd) | lcapy.C(Cp) | lcapy.R(Rpolar)) + lcapy.C(Cc)) |
#            lcapy.C(Chemt))

# circuit.draw()
R_noise = circuit.noisy()
R_noise.draw()

in_noise = 1e-9*frequency
circuit = circuit.subs({'Cd':Cd})
imped = circuit[2].Z(lcapy.f).evaluate(frequency)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(frequency, abs(imped), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)

# vc = circuit.R.v.evaluate(frequency)