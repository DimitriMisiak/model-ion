#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:59:15 2020

@author: filippini

General Hemt model with a single electrode, with a arbitrary impedance
and voltage noise generator (thevenin noise circuit) on the input node E.
The Hemt has a resistive feedback.
"""

import lcapy
from lcapy import Circuit, s, f, j, t, omega, expr, oo, symbol
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

cir = Circuit("""
    O E 0; down
    W 0 0_g; down=0.2, cground                   
    
    W E 1_r; right
    Rpolar 1_r 0_r; down
    Vj 0_r 0_j s {Vj}; down, l={e_J,polar}, color=red
    Cc 1_r 1_fb; right=1.5
    Cfb 1_fb 0_fb; down
    W 1_fb 1_h; right
    Chemt 1_h 0_h; down
    W 1_h 1_s; right

    In 1_s 0_s s {In}; down, color=red, l=i_n
    
    W 0 0_j; right
    W 0_j 0_fb; right
    W 0_fb 0_h; right
    W 0_h 0_s; right
    

    ; draw_nodes=primary
    ; label_nodes=alpha
    ; style=european
    ; label_ids=False
""")

cir.draw()

source_list = ['In', 'Vj']

Z = cir.impedance('E')

E_dict = dict()
for source in source_list:
    cir_aux = cir.kill_except(source)
    E_thevenin = cir_aux['E'].V(s)
    E_dict[source] = E_thevenin
    
I_dict = {k:v/Z for k,v in E_dict.items()}
I_norm_dict ={k:v/k for k,v in I_dict.items()}
