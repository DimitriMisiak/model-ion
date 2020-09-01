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


# ### electron-holes pairs on 2
# cir = Circuit("""
#     W 0 0_0; down=0.2, cground          
    
#     C11 1 0; down
#     C22 2 0_2; down
#     C12 2 1; right, size=2
    
#     Iq 2 1 {Npair*qe*(u(t) - u(t-1))}; right, offset=1
    
#     W 0_2 0; right
#     ; draw_nodes=connections
#     ; label_nodes=none
#     ; style=european
# """)

# cir = Circuit("""
              

#     C12 A B; rotate=-45, size=2
#     C23 B C; rotate=-135, size=2
#     C34 C D; rotate=135, size=2
#     C14 D A; rotate=45, size=2

#     C24 D 1_0; right
#     W 1_0 B; right
#     C13 A 1_1; down
#     W 1_1 C; down


#     C11 0_1 A; right
#     W 0_1 0_4; down
#     C44 0_4 D; right
#     W 0_4 0_3; down
#     C33 0_3 C; right
#     W 0_3 0_2; down
    
#     W 0_2 0; down
#     C22 0_2 2_0; right
#     W B 2_0; down
    
#     W 0 0_0; down=0.2, cground

#     ; draw_nodes=connections
#     ; label_nodes=all
#     ; style=european
# """)

cir = Circuit("""
              

    Z12 A B; rotate=-45, size=2
    Z23 B Z; rotate=-135, size=2
    Z34 Z D; rotate=135, size=2
    Z14 D A; rotate=45, size=2

    Z24 D 4_2; right
    W 4_2 B; right
    Z13 A 1_3; down
    W 1_3 Z; down

    W 1 A; down
    W 1_1 1; right
    W 4 D; right
    W Z 3; down
    W B 2; right

    Z33 3 0; down
    Z44 4 0_4; down
    Z22 2 0_2; down
    Z11 1_1 0_1; down
    
    W 0_1 0_4; right
    W 0_4 0; right
    W 0 0_2; right
    W 0 0_0; down=0.2, cground
    
    W 2 2_2; right

    Iq 2 0_2 {Npair*qe*(u(t) - u(t-1))}; down, offset=1

    W 0_2 2_3; right
    
    ; draw_nodes=connections
    ; label_nodes=alpha
    ; style=european
""")


cir.draw()

cir_simple = cir.kill()

cir_eval = cir.subs(mutual_dict)

# # pulse_dict = {
# #     'Iq': expr('e*(u(t) - u(t-1))')
# # }

# cir_eval = cir.subs(mutual_dict)

# t_array = np.linspace(-10, 10, 1000)

# # # pulse
# # cir_eval.Iq.I(t).plot(t_array)

# # cir_eval[1].V(t).plot(t_array)

# V1 = (cir_eval[1].V(t))(2).evalf()
# V2 = (cir_eval[2].V(t))(2).evalf()
