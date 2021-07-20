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
from tqdm import tqdm
import sympy as sy

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

Zhemt_Cfb = "((Cc + Cfb + Chemt)/(Cc*(Cfb + Chemt)))/(s + (Cc + Cfb + Chemt)/(Cc*Cfb*Rpolar + Cc*Chemt*Rpolar))"

In_Cfb = "Cc*In/(Cc + Cfb + Chemt)"
Ij_Cfb = "Vj/Rpolar"

Vn = '(e0**2 +(ea/f**0.5)**2 + (eb/f)**2)**0.5 '
In = '(i0**2 +(ia*(f)**0.5)**2 + (ib*f)**2)**0.5 '
Vj = '(4*k*T*Rpolar)**0.5 '
Vfb = '(4*k*T*Rfb)**0.5 '

InN = f'Cc*{In}/(Cc + Cfb + Chemt)'
IjN = f'{Vj}/Rpolar'

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

# =============================================================================
# NETLIST
# =============================================================================
cir = Circuit(fr"""

    W 0 GND; down=0.2, cground,  
    
    W 3_0 3_ej0; right,  
    W 3_ej0 3_i0; right,  
    W 3_i0 3_r0; right,  
    W 3_r0 3_p0; right,  
    W 3_p0 3_d0; right,  
    W 3_d0 0; right=3,  

    W 0 4_d0; right=3,  
    W 4_d0 4_p0; right,  
    W 4_p0 4_r0; right,  
    W 4_r0 4_i0; right,  
    W 4_i0 4_ej0; right,  
    W 4_ej0 4_0; right,  

    W 1_0 1_ej0; right,  
    W 1_ej0 1_i0; right,  
    W 1_i0 1_r0; right,  
    W 1_r0 1_p0; right,  
    W 1_p0 1_d0; right,  

    W 2_d0 2_p0; right,  
    W 2_p0 2_j0; right,  
    W 2_j0 2_fb0; right,  
    W 2_fb0 2_h0; right,  
    W 2_h0 2_s0; right,  
    W 2_s0 2_0; right,  
    
    W 1_0 3_0; down=3,  
    W 2_0 4_0; down=3,  

    Ij1 1_ej 1_ej0 noise {{{IjN}}}; down, l=i_{{j,N}}^A, color=red
    W 1_ej 1_i; right
    In1 1_i 1_i0 noise {{{InN}}}; down, l=i_{{n,N}}^A, color=red
    W 1_i 1_r
    Zhemt1 1_r 1_r0 {{{Zhemt_Cfb}}}; down, l=Z_{{hemt}}
    W 1_r 1_p; right
    Cp1 1_p 1_p0 Cp; down, l=C_{{p}}
    W 1_p 1_d; right
    C11 1_d 1_d0; down, l=C_{{11}}^m
    W 1_d 1_A; right=1
    W 1_A A; right=1
    
    W B 2_B; right=1
    W 2_B 2_d; right=1
    C22 2_d 2_d0; down, l=C_{{22}}^m
    W 2_d 2_p; right
    Cp2 2_p 2_p0 Cp; down, l=C_{{p}}
    W 2_p 2_r; right
    Rpolar2 2_r 2_r0 Rpolar; down, l=R_{{ polar.}}
    Vj2 2_r0 2_j0 noise {{{Vj}}}; down, l=e_J^B, color=red
    Cc2 2_r 2_fb Cc; right=1.5, l=C_{{c}}
    Cfb2 2_fb 2_fb0 Cfb; down, l=C_{{fb}}
    W 2_fb 2_h; right
    Chemt2 2_h 2_h0 Chemt; down, l=C_{{hemt}}
    Vn2 2_h 2_s noise {{{Vn}}}; right, color=red, l=e_n^B
    In2 2_s 2_s0 noise {{{In}}}; down, color=red, l=i_n^B
    W 2_s S_B; right, endarrow=otri, draw_nodes=none    
    
    W D 4_D; right=1
    W 4_D 4_d; right=1
    C44 4_d 4_d0; down, l=C_{{44}}^m
    W 4_d 4_p; right
    Cp4 4_p 4_p0 Cp; down, l=C_{{p}}
    W 4_p 4_r; right
    Zhemt4 4_r 4_r0 {{{Zhemt_Cfb}}}; down, l=Z_{{hemt}}
    W 4_r 4_i; right
    In4 4_i 4_i0 noise {{{InN}}}; down, l=i_{{n,N}}^C, color=red 
    W 4_i 4_ej; right
    Ij4 4_ej 4_ej0 noise {{{IjN}}}; down, l=i_{{j,N}}^C, color=red

    Ij3 3_ej 3_ej0 noise {{{IjN}}}; down, l=i_{{j,N}}^C, color=red
    W 3_ej 3_i; right
    In3 3_i 3_i0 noise {{{InN}}}; down, l=i_{{n,N}}^C, color=red
    W 3_i 3_r; right
    Zhemt3 3_r 3_r0 {{{Zhemt_Cfb}}}; down, l=Z_{{hemt}}
    W 3_r 3_p; right
    Cp3 3_p 3_p0 Cp; down, l=C_{{p}}
    W 3_p 3_d; right
    C33 3_d 3_d0; down, l=C_{{33}}^m
    W 3_d 3_C; right=1
    W 3_C C; right=1

    C13 A C; down=4, l=C_{{13}}^m
    C12 A B; right=4, l=C_{{12}}^m
    C24 B D; down=4, l=C_{{24}}^m
    C34 C D; right=4, l=C_{{34}}^m
    
    W A 14; rotate=-45
    C14 14 D; free, scale=0.75, l=C_{{14}}^m

    W B 23; rotate=-135
    C23 23 C; free, scale=0.75, l=C_{{23}}^m
    
    Iequator 1_A 3_C {{Npair*qe*delta(t)}}; down, color=blue, l_=I_s^{{equator}}
    Iveto 1_A 2_B {{Npair*qe*delta(t)}}; right, color=blue, offset=1, l_=I_s^{{veto,top}}
    Ibulk 2_B 4_D {{Npair*qe*delta(t)}}; down, color=blue, l^=I_s^{{bulk}}
        
    ; draw_nodes=connections
    ; label_nodes=alpha
    ; style=european
    ; help_lines=0

""")    

cir.draw()

print(cir.netlist())