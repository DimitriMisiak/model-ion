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
import lcapy
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

maxwell_matrix = np.loadtxt('fid38_maxwell.txt', comments='%')

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
    'C13': mutual_matrix[0,2], #F
    'C14': mutual_matrix[0,3], #F
    'C22': mutual_matrix[1,1], #F
    'C23': mutual_matrix[1,2], #F
    'C24': mutual_matrix[1,3], #F
    'C33': mutual_matrix[2,2], #F
    'C34': mutual_matrix[2,3], #F
    'C44': mutual_matrix[3,3], #F
}

eval_dict = {
    'Cd': 20e-12, #F 
    'Cp': 10e-12, #F
    'Rpolar': 1e10, #Ohms
    'Cc': 2e-9, #F
    'Chemt': 5e-12, #F
    'Rfb': 1e10, #1e10, #F
    'Cfb': 3e-12, #F
    'qe': 1.60e-19, #C,
    'Npair': 1e3/3, # number of pairs in 1keV
    'dt': 1e-5
}
eval_dict.update(mutual_dict)

cir = Circuit(r"""

    W 0 GND; down=0.2, cground, ultra thick
    
    W 3_0 3_h0; right=2, ultra thick
    W 3_h0 3_fb0; right, ultra thick
    W 3_fb0 3_p0; right, ultra thick
    W 3_p0 3_d0; right, ultra thick
    W 3_d0 0; right=2, ultra thick

    W 0 4_d0; right=2, ultra thick
    W 4_d0 4_p0; right, ultra thick
    W 4_p0 4_fb0; right, ultra thick
    W 4_fb0 4_h0; right, ultra thick
    W 4_h0 4_0; right=2, ultra thick

    W 1_0 1_h0; right=2, ultra thick
    W 1_h0 1_fb0; right, ultra thick
    W 1_fb0 1_p0; right, ultra thick
    W 1_p0 1_d0; right, ultra thick

    W 2_d0 2_p0; right, ultra thick
    W 2_p0 2_fb0; right, ultra thick
    W 2_fb0 2_h0; right, ultra thick
    W 2_h0 2_0; right=2, ultra thick
    
    W 1_0 3_0; down=2, ultra thick
    W 2_0 4_0; down=2, ultra thick

    W S_A 1_h; right, startarrow=otri, draw_nodes=none
    Chemt1 1_h 1_h0 Chemt; down, l=C_{hemt}
    W 1_h 1_fb; right
    Cfb1 1_fb 1_fb0 Cfb; down, l=C_{fb}
    Cc1 1_fb 1_p Cc; right=1.5, l=C_{c}
    Cp1 1_p 1_p0 Cp; down, l=C_{p}
    W 1_p 1_d; right
    C11 1_d 1_d0; down, l=C_{11}
    W 1_d A; right
    
    W B 2_d; right
    C22 2_d 2_d0; down, l=C_{22}
    W 2_d 2_p; right
    Cp2 2_p 2_p0 Cp; down, l=C_{p}
    Cc2 2_p 2_fb Cc; right=1.5, l=C_{c}
    Cfb2 2_fb 2_fb0 Cfb; down, l=C_{fb}
    W 2_fb 2_h; right
    Chemt2 2_h 2_h0 Chemt; down, l=C_{hemt}
    W 2_h S_B; right, endarrow=otri, draw_nodes=none    
    
    W D 4_d; right
    C44 4_d 4_d0; down, l=C_{44}
    W 4_d 4_p; right
    Cp4 4_p 4_p0 Cp; down, l=C_{p}
    Cc4 4_p 4_fb Cc; right=1.5, l=C_{c}
    Cfb4 4_fb 4_fb0 Cfb; down, l=C_{fb}
    W 4_fb 4_h; right
    Chemt4 4_h 4_h0 Chemt; down, l=C_{hemt}
    W 4_h S_D; right, endarrow=otri, draw_nodes=none    
    
    W S_C 3_h; right, startarrow=otri, draw_nodes=none
    Chemt3 3_h 3_h0 Chemt; down, l=C_{hemt}
    W 3_h 3_fb; right
    Cfb3 3_fb 3_fb0 Cfb; down, l=C_{fb}
    Cc3 3_fb 3_p Cc; right=1.5, l=C_{c}
    Cp3 3_p 3_p0 Cp; down, l=C_{p}
    W 3_p 3_d; right
    C33 3_d 3_d0; down, l=C_{33}
    W 3_d C; right

    C13 A C; down=2
    C12 A B; right=2
    C24 B D; down=2
    C34 C D; right=2
    
    W A 14; rotate=-45
    C14 14 D; free, scale=0.75

    W B 23; rotate=-135
    C23 23 C; free, scale=0.75
    
    ; draw_nodes=connections
    ; label_nodes=alpha
    ; style=european
    ; help_lines=1
    
    ;;\node[blue,draw,dashed,inner sep=10mm, fit=(1_d) (2_d) (0), label=Detector] {};
    ;;\node[blue,draw,dashed,inner sep=5mm, fit=(S_A) (1_h0) (Cc1), label=Hemt 1] {};
    ;;\node[blue,draw,dashed,inner sep=5mm, fit=(S_B) (2_h0) (Cc2), label=Hemt 2] {};
    ;;\node[blue,draw,dashed,inner sep=5mm, fit=(S_C) (3_h0) (Cc3), label=Hemt 3] {};
    ;;\node[blue,draw,dashed,inner sep=5mm, fit=(S_D) (4_h0) (Cc4), label=Hemt 4] {};
""")


e_n = lcapy.expr('(e0**2 +(ea/f**0.5)**2 + (eb/f)**2)**0.5 ')
i_n = lcapy.expr('(i0**2 +(ia*(f)**0.5)**2 + (ib*f)**2)**0.5 ')
e_j = lcapy.expr('(4*k*T*Rpolar)**0.5 ')
e_fb = lcapy.expr('(4*k*T*Rfb)**0.5 ')

noise_source_dict = {
    'Vj': e_j,
    'Vfb': e_fb,
    'Vdaq': 0,
    'In': i_n,
    'Vn': e_n,
}

noise_eval_dict = {
    'T': 20e-3, #K
    'k': 1.38e-23, #J/K,
    'e0': 0.22e-9, #V
    'ea': 36e-9, #V*Hz**0.5
    'eb': 0, # V*Hz
    'i0': 4.0e-23, #A
    'ia': 2.6e-18, #A/Hz**0.5
    'ib': 0, #A/Hz
}
noise_eval_dict.update(eval_dict)

cir.draw()

cir_eval =cir.subs(eval_dict)