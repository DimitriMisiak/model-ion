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

It checks out with the calculation on paper :D
"""

from lcapy import Circuit, s, f, j, omega
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

an_dict = {
    'R1': 1e10, #Omhs
    'R2': 1e10, #Ohms
    'C11': 5e-12, #F
    'C12': 5e-12, #F
    'C22': 5e-12, #F
}

cir = Circuit("""
    W 0 0_0; down, ground, size=0.1           
    
    C11 1 0; down
    C22 2 0_2; down
    C12 2 1; right, size=2
    
    W 1 1_1; right
    R1 1_1 0_11; down
    W 2_1 2; right
    R2 2_1 0_21; down
    
    Vj1 0_11 0_111 {Vj1(s)}; down
    Vj2 0_21 0_211 {Vj2(s)}; down

    W 0 0_111; right
    W 0_2 0; right
    W 0_211 0_2; right
    ; draw_nodes=connections
""")

cir.draw()

cir2 = cir.kill_except('Vj2')
e2 = cir2.Vj2.V(s)
s2 = cir2[1].V(s)
H2 = (s2/e2).simplify()

cir1 = cir.kill_except('Vj1')
e1 = cir1.Vj1.V(s)
s1 = cir1[1].V(s)
H1 = (s1/e1).simplify()

A = (H2/H1).frequency_response()

A_eval = A.subs(an_dict)

A_abs = A_eval.abs
A_phase = A_eval.phase

f_array = 10**np.linspace(-1, 2, 1000)
abs_array = A_abs.evaluate(f_array)
phase_array = A_phase.evaluate(f_array)

fig, ax = plt.subplots(nrows=2, sharex=True)
ax[0].set_title("""
    Ratio of the response in node 1 of
    the Johnson noise of R2 over the Johnson noise of R1
""")

ax[0].plot(f_array, abs_array)

ax[0].set_ylabel('Absolute Value [Fraction]')
ax[0].set_xlabel('Frequency [Hz]')

ax[1].plot(f_array, phase_array)
ax[1].set_xscale('log')

ax[1].set_ylabel('Phase [Radians]')
ax[1].set_xlabel('Frequency [Hz]')

for a in ax:
    a.grid()

fig.tight_layout()
