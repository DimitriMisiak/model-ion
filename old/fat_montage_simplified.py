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

Zhemt_Cfb = "((Cc + Cfb + Chemt)/(Cc*(Cfb + Chemt)))/(s + (Cc + Cfb + Chemt)/(Cc*Cfb*Rpolar + Cc*Chemt*Rpolar))"

In_Cfb = "Cc*In/(Cc + Cfb + Chemt)"
Ij_Cfb = "Vj/Rpolar"

cir = Circuit(fr"""

    W 0 GND; down=0.2, cground, ultra thick
    
    W 3_0 3_ej0; right, ultra thick
    W 3_ej0 3_i0; right, ultra thick
    W 3_i0 3_r0; right, ultra thick
    W 3_r0 3_p0; right, ultra thick
    W 3_p0 3_d0; right, ultra thick
    W 3_d0 0; right=3, ultra thick

    W 0 4_d0; right=3, ultra thick
    W 4_d0 4_p0; right, ultra thick
    W 4_p0 4_r0; right, ultra thick
    W 4_r0 4_i0; right, ultra thick
    W 4_i0 4_ej0; right, ultra thick
    W 4_ej0 4_0; right, ultra thick

    W 1_0 1_ej0; right, ultra thick
    W 1_ej0 1_i0; right, ultra thick
    W 1_i0 1_r0; right, ultra thick
    W 1_r0 1_p0; right, ultra thick
    W 1_p0 1_d0; right, ultra thick

    W 2_d0 2_p0; right, ultra thick
    W 2_p0 2_j0; right, ultra thick
    W 2_j0 2_fb0; right, ultra thick
    W 2_fb0 2_h0; right, ultra thick
    W 2_h0 2_s0; right, ultra thick
    W 2_s0 2_0; right, ultra thick
    
    W 1_0 3_0; down=3, ultra thick
    W 2_0 4_0; down=3, ultra thick

    Ij1 1_ej 1_ej0 noise; down, l=i_{{j,N}}, color=red
    W 1_ej 1_i; right
    In1 1_i 1_i0 noise; down, l=i_{{n,N}}, color=red
    W 1_i 1_r
    Zhemt1 1_r 1_r0 {{{Zhemt_Cfb}}}; down, l=Z_{{hemt}}
    W 1_r 1_p; right
    Cp1 1_p 1_p0 Cp; down, l=C_{{p}}
    W 1_p 1_d; right
    C11 1_d 1_d0; down, l=C_{{11}}
    W 1_d 1_A; right=1
    W 1_A A; right=1
    
    W B 2_B; right=1
    W 2_B 2_d; right=1
    C22 2_d 2_d0; down, l=C_{{22}}
    W 2_d 2_p; right
    Cp2 2_p 2_p0 Cp; down, l=C_{{p}}
    W 2_p 2_r; right
    Rpolar2 2_r 2_r0 Rpolar; down, l=R_{{load}}
    Vj2 2_r0 2_j0 noise ej; down, l=e_J, color=red
    Cc2 2_r 2_fb Cc; right=1.5, l=C_{{c}}
    Cfb2 2_fb 2_fb0 Cfb; down, l=C_{{fb}}
    W 2_fb 2_h; right
    Chemt2 2_h 2_h0 Chemt; down, l=C_{{hemt}}
    Vn2 2_h 2_s noise Vn; right, color=red, l=e_n
    In2 2_s 2_s0 noise In; down, color=red, l=i_n
    W 2_s S_B; right, endarrow=otri, draw_nodes=none    
    
    W D 4_D; right=1
    W 4_D 4_d; right=1
    C44 4_d 4_d0; down, l=C_{{44}}
    W 4_d 4_p; right
    Cp4 4_p 4_p0 Cp; down, l=C_{{p}}
    W 4_p 4_r; right
    Zhemt4 4_r 4_r0 {{{Zhemt_Cfb}}}; down, l=Z_{{hemt}}
    W 4_r 4_i; right
    In4 4_i 4_i0 noise; down, l=i_{{n,N}}, color=red 
    W 4_i 4_ej; right
    Ij4 4_ej 4_ej0 noise; down, l=i_{{j,N}}, color=red

    Ij3 3_ej 3_ej0 noise; down, l=i_{{j,N}}, color=red
    W 3_ej 3_i; right
    In3 3_i 3_i0 noise; down, l=i_{{n,N}}, color=red
    W 3_i 3_r; right
    Zhemt3 3_r 3_r0 {{{Zhemt_Cfb}}}; down, l=Z_{{hemt}}
    W 3_r 3_p; right
    Cp3 3_p 3_p0 Cp; down, l=C_{{p}}
    W 3_p 3_d; right
    C33 3_d 3_d0; down, l=C_{{33}}
    W 3_d 3_C; right=1
    W 3_C C; right=1

    C13 A C; down=4
    C12 A B; right=4
    C24 B D; down=4
    C34 C D; right=4
    
    W A 14; rotate=-45
    C14 14 D; free, scale=0.75

    W B 23; rotate=-135
    C23 23 C; free, scale=0.75
    
    Iequator 1_A 3_C {{Npair*qe*delta(t)}}; down, color=blue, l=I_{{equator}}
    Iveto 1_A 2_B {{Npair*qe*delta(t)}}; right, color=blue, offset=1, l^=I_{{veto}}
    Ibulk 2_B 4_D {{Npair*qe*delta(t)}}; down, color=blue, l^=I_{{bulk}}
        
    ; draw_nodes=connections
    ; label_nodes=alpha
    ; style=european
    ; help_lines=1

""")    

e_n = lcapy.expr('(e0**2 +(ea/f**0.5)**2 + (eb/f)**2)**0.5 ')
i_n = lcapy.expr('(i0**2 +(ia*(f)**0.5)**2 + (ib*f)**2)**0.5 ')
e_j = lcapy.expr('(4*k*T*Zhemt)**0.5 ')
e_fb = lcapy.expr('(4*k*T*Rfb)**0.5 ')
i_n_Cfb = lcapy.expr(In_Cfb)
i_j_Cfb = lcapy.expr(Ij_Cfb)

noise_source_dict = {
    'Vj': e_j,
    'Vfb': e_fb,
    'Vdaq': 0,
    'In': i_n,
    'Vn': e_n,
    'In1': i_n_Cfb,
    'In3': i_n_Cfb,
    'In4': i_n_Cfb,
    'Ij1': i_j_Cfb,
    'Ij3': i_j_Cfb,
    'Ij4': i_j_Cfb,
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

cir_eval = cir.subs(eval_dict)

# Z = cir.Zhemt1.Z.symbols['Zhemt']
# ZZ = lcapy.expr("25100000000000/(953*(s + 2510/953))")

cir_kill = cir_eval.kill()
cir_veto = cir_eval.kill_except('Iveto')
cir_bulk = cir_eval.kill_except('Ibulk')
cir_equator = cir_eval.kill_except('Iequator')

cir_list = [cir_veto, cir_bulk, cir_equator]


# A = cir_bulk['S_B'].V(s)
# A_eval = A.subs({
#     Z:25100000000000/(953*(s + 2510/953))
# }).evalf()
# print('Done!\n')
# print(A)

#%%
# # =============================================================================
# # SIGNAL GENERATION
# # =============================================================================
# cir_signal = cir_bulk.subs(eval_dict)
# freq_array = 10**np.linspace(0, 4, int(1e3))

# # frequency plot
# SIGNAL = abs(cir_signal['S_B'].V(j*omega)(f).evaluate(freq_array))

# fig, ax = plt.subplots(nrows=2, sharex=True)
# ax[1].loglog(freq_array, SIGNAL, label='Output, Node S')
# ax[1].set_ylabel('Voltage [V]')
# ax[1].set_xlabel('Frequency [Hz]')

# for a in ax:
#     a.grid()
#     a.legend()

#%%
from tqdm import tqdm
# =============================================================================
# NOISE PROPAGATION
# =============================================================================
noise_sources = [source for source in cir.sources if cir[source].is_noisy]

cir_noise = cir_eval.kill_except(*noise_sources)

tot_lab = 'Tot'
noise_expr_dict = dict()
for source in tqdm(noise_sources) :  #+ [tot_lab,]
    
    if source == tot_lab:
        cir_source = cir_noise
    else:
        cir_source = cir_noise.kill_except(source)

    voutput_source = cir_source['S_B'].V.n
    
    noise_expr_dict[source] = voutput_source


fig, ax = plt.subplots()

for source, noise_expr in noise_expr_dict.items():

    noise_funk = noise_expr.subs(noise_source_dict).subs(noise_eval_dict)
    # hack
    noise_fun = noise_funk(omega)(f) 
    noise_array = noise_fun.evaluate(freq_array)
    
    if source == tot_lab:
        NOISE = noise_array
        ax.loglog(freq_array, noise_array, label=source, lw=2, color='k')
    else:
        ax.loglog(freq_array, noise_array, label=source)

ax.set_ylabel(r'LPSD [V/$\sqrt{Hz}$]')
ax.set_xlabel(r'Frequency [Hz]')
ax.legend()
ax.grid()