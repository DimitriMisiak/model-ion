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

from lcapy import Circuit, s, f, j, t, omega, expr, oo, symbol
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

maxwell_matrix = np.loadtxt('planar38_maxwell.csv', comments='%', delimiter=',')

def maxwell_to_mutual(mat):
    """
    From the Maxwell capacitance matrix, return the mutual capacitance matrix.
    From the Mutual capacitance matrix, return the Maxwell capacitance matrix.
    The transformation is its own inverse, that is why.
    """
    ret = np.zeros(mat.shape)
    
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if i == j:
                ret[i,j] = sum(mat[i])
            else:
                ret[i,j] = -mat[i,j]
                
    return ret

mutual_matrix = maxwell_to_mutual(maxwell_matrix)

mutual_dict = {
    'C11': mutual_matrix[0,0], #F
    'C12': mutual_matrix[0,1], #F
    'C22': mutual_matrix[1,1], #F
    'qe': 1.60e-19, #C,
    'Npair': 1e3/3, # number of pairs in 1keV
}

qe = symbol('qe', real=True)
Npair = symbol('Npair', real=True)

### holes only on 1
# cir = Circuit("""
#     W 0 0_0; down=0.2, cground          
    
#     C11 1 0; down
#     C22 2 0_2; down
#     C12 2 1; right, size=2
    
#     Iq 1 0 {Npair*qe*(u(t) - u(t-1))}; down, offset=1
    
#     W 0_2 0; right
#     ; draw_nodes=connections
#     ; label_nodes=none
#     ; style=european
# """)


### electron-holes pairs on 2
cir = Circuit("""
    W 0 0_0; down=0.2, cground          
    
    C11 1 0; down
    C22 2 0_2; down
    C12 2 1; right, size=2
    
    Iq 2 1 {Npair*qe*(u(t) - u(t-1))}; right, offset=1
    
    W 0_2 0; right
    ; draw_nodes=connections
    ; label_nodes=none
    ; style=european
""")


cir.draw()

# pulse_dict = {
#     'Iq': expr('e*(u(t) - u(t-1))')
# }

cir_eval = cir.subs(mutual_dict)

t_array = np.linspace(-10, 10, 1000)

# # pulse
# cir_eval.Iq.I(t).plot(t_array)

# cir_eval[1].V(t).plot(t_array)

V1 = (cir_eval[1].V(t))(2).evalf()
V2 = (cir_eval[2].V(t))(2).evalf()

T = cir.transient()

K = cir.kill()

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
