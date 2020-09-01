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
    'qe': 1.60e-19, #C,
    'Npair': 1e3/3, # number of pairs in 1keV
}

qe = symbol('qe', real=True)
Npair = symbol('Npair', real=True)

A = symbol('A', function=True)


funk_dict = {
    A: 'Npair*qe*(u(t) - u(t-1))'
}

cir = Circuit("""
    C12 A B; rotate=-45, size=2
    C23 B C; rotate=-135, size=2
    C34 C D; rotate=135, size=2
    C14 D A; rotate=45, size=2

    C24 D 4_2; right
    W 4_2 B; right
    C13 A 1_3; down
    W 1_3 C; down

    W 1 A; down
    W 1_1 1; right
    W 4 D; right
    W C 3; down
    W B 2; right

    C33 3 0; down
    C44 4 0_4; down
    C22 2 0_2; down
    C11 1_1 0_1; down
    
    W 0_1 0_4; right
    W 0_4 0; right
    W 0 0_2; right
    W 0 0_0; down=0.2, cground
    
    W 2 2_2; right=2

    IB 2 0_2 {Npair*qe*(u(t) - u(t-1))}; down, offset=1, color=blue, l^=I_{B}



    IA 1_1 0_1 {Npair*qe*(u(t) - u(t-1))}; down, offset=-1, color=blue, l_=I_{A}

    IAB A 1_2 {Npair*qe*(u(t) - u(t-1))}; scale=0.5, right, color=blue, l^=I_{veto}
    W 1_2 B; down, color=blue

    W D 4_B; down=2, color=blue
    W 4_B 4_BB; right, color=blue
    IBD 4_BB 2_D {Npair*qe*(u(t) - u(t-1))}; scale=0.5, right, color=blue, l_=I_{bulk}
    W B 2_D; down, color=blue

    W 1_C A; right=2, color=blue
    IAC 1_C 3_A {Npair*qe*(u(t) - u(t-1))}; scale=0.5, down, color=blue, l_=I_{equator}
    W 3_A 3_AA; down, color=blue
    W 3_AA C; right, color=blue


    W 0_2 2_3; right
    
    ; draw_nodes=connections
    ; label_nodes=alpha
    ; style=european
""")

# Rp 2 0_2 ; down, offset=2

cir.draw()

cir_kill = cir.kill()

cir_eval = cir.subs(mutual_dict)

cir_ib = cir_eval.kill_except('IB')
cir_ia = cir_eval.kill_except('IA')
cir_veto = cir_eval.kill_except('IAB')
cir_bulk = cir_eval.kill_except('IBD')
cir_equator = cir_eval.kill_except('IAC')

cir_list = [cir_ib, cir_ia, cir_veto, cir_bulk, cir_equator]

# ### Signal Generation

# volt_list = list()
# for ctt in cir_list:
#     v_list = list()
#     for i in range(1,5):
#         v = (ctt[i].V(t))(2).evalf()
#         v_list.append(v)
    
#     volt_list.append(v_list)

# v_array = np.array(volt_list).astype(float)

# print(v_array*1e6) #µV/keV