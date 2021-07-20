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

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1

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

alpha = 1

res_alpha = list()

# alpha_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2]
capa_list = np.array([0, 2.5, 5, 7.5, 10, 12.25, 15, 17.25, 20, 25, 30]) * 1e-12 #F

for capa_cabling in tqdm(capa_list):
    
    mutual_matrix = maxwell_to_mutual(maxwell_matrix) * alpha
    
    beta = 1
    mutual_dict = {
        'C11': mutual_matrix[0,0], #F
        'C12': mutual_matrix[0,1]*beta, #F
        'C22': mutual_matrix[1,1], #F
    }

    eval_dict = {
        'Cp': capa_cabling, #F
        'Rpolar': 1e10, #Ohms
        'Cc': 2e-9, #F
        'Chemt': 5e-12, #F
        'Rfb': 1e10, #1e10, #F
        'Cfb': 3e-12, #F
        'qe': 1.60e-19, #C,
        'Npair': 1e3/3, # number of pairs in 1keV
        'dt': 1e-5
    }
    
    # eval_dict = {
    #     'Cp': 5e-12, #F
    #     'Rpolar': 1e10, #Ohms
    #     'Cc': 2e-9, #F
    #     'Chemt': 5e-12, #F
    #     'Rfb': 1e10, #1e10, #F
    #     'Cfb': 3e-12, #F
    #     'qe': 1.60e-19, #C,
    #     'Npair': 1e3/3, # number of pairs in 1keV
    #     'dt': 1e-5
    # }
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
    
    # noise_eval_dict = {
    #     'T': 20e-3, #K
    #     'k': 1.38e-23, #J/K,
    #     'e0': 0.22e-9, #V
    #     'ea': 36e-9, #V*Hz**0.5
    #     'eb': 0, # V*Hz
    #     'i0': 4.0e-23, #A
    #     'ia': 2.6e-18, #A/Hz**0.5
    #     'ib': 0, #A/Hz
    # }
    
    noise_eval_dict = {
        'T': 10e-3, #K
        'k': 1.38e-23, #J/K,
        'e0': 2.56e-10, #V
        'ea': 2.965e-8, #V*Hz**0.5
        'eb': 5.343e-9, # V*Hz
        'i0': 4.0e-23, #A
        'ia': 2.6e-18, #A/Hz**0.5
        'ib': 0, #A/Hz
    }
    
    noise_eval_dict.update(eval_dict)
    
    # =============================================================================
    # NETLIST
    # =============================================================================
    cir = Circuit(fr"""
    
        W 0 GND; down=0.2, cground, ultra thick
    
        W 1_s0 1_h0; right, ultra thick
        W 1_h0 1_fb0; right, ultra thick
        W 1_fb0 1_j0; right, ultra thick
        W 1_j0 1_p0; right, ultra thick
        W 1_p0 1_d0; right, ultra thick
        W 1_d0 0; right=2, ultra thick
    
        W 0 2_d0; right=2, ultra thick
        W 2_d0 2_p0; right, ultra thick
        W 2_p0 2_j0; right, ultra thick
        W 2_j0 2_fb0; right, ultra thick
        W 2_fb0 2_h0; right, ultra thick
        W 2_h0 2_s0; right, ultra thick
    
        W S_A 1_s; right, startarrow=otri, draw_nodes=none
        In1 1_s 1_s0 noise {{{In}}}; down, color=red, l=i_n
        Vn1 1_s 1_h noise {{{Vn}}}; right, color=red, l=e_n
        Chemt1 1_h 1_h0 Chemt; down, l=C_{{hemt}}
        W 1_h 1_fb; right
        Cfb1 1_fb 1_fb0 Cfb; down, l=C_{{fb}}
        Cc1 1_fb 1_r Cc; right=1.5, l=C_{{c}}
        Rpolar1 1_r 1_r0 Rpolar; down, l=R_{{load}}
        Vj1 1_r0 1_j0 noise {{{Vj}}}; down, l=e_J, color=red
        W 1_r 1_p; right
        Cp1 1_p 1_p0 Cp; down, l=C_{{p}}
        W 1_p 1_d; right
        C11 1_d 1_d0; down, l=C_{{11}}
        W 1_d 1_A; right
        W 1_A A; right
        
        W B 2_B; right=1
        W 2_B 2_d; right=1
        C22 2_d 2_d0; down, l=C_{{22}}
        W 2_d 2_p; right
        Cp2 2_p 2_p0 Cp; down, l=C_{{p}}
        W 2_p 2_r; right
        Rpolar2 2_r 2_r0 Rpolar; down, l=R_{{load}}
        Vj2 2_r0 2_j0 noise {{{Vj}}}; down, l=e_J, color=red
        Cc2 2_r 2_fb Cc; right=1.5, l=C_{{c}}
        Cfb2 2_fb 2_fb0 Cfb; down, l=C_{{fb}}
        W 2_fb 2_h; right
        Chemt2 2_h 2_h0 Chemt; down, l=C_{{hemt}}
        Vn2 2_h 2_s noise {{{Vn}}}; right, color=red, l=e_n
        In2 2_s 2_s0 noise {{{In}}}; down, color=red, l=i_n
        W 2_s S_B; right, endarrow=otri, draw_nodes=none    
    
        C12 A B; right=4
        
        Ibulk 1_A 2_B {{Npair*qe*delta(t)}}; right, color=blue, offset=1, l^=I_{{veto}}
            
        ; draw_nodes=connections
        ; label_nodes=alpha
        ; style=european
        ; help_lines=0
        ; font=\normalsize
    
    """)    
    
    cir.draw()
    
    ref_point = 'S_B'
    
    signal_sources = [source for source in cir.sources if not cir[source].is_noisy]
    noise_sources = [source for source in cir.sources if cir[source].is_noisy]
    
    freq_array = 10**np.linspace(0, 4, int(1e3))
    
    #%%
    # =============================================================================
    # STEP GENERATION Generation
    # =============================================================================
    cir_step = cir.subs({'Rpolar':oo}).subs(eval_dict)
    
    step_dict = dict()
    
    for source in signal_sources:
        cir_aux = cir_step.kill_except(source)
        step_dict[source] = (cir_aux[ref_point].V(s)*s).evalf()
    
    Y_dict = dict()
    Y_dict['B'] = (cir_step.admittance('B')/s).evalf()
    Y_dict['S_B'] = (cir_step.admittance('S_B')/s).evalf()
    
    #%%
    # =============================================================================
    # SIGNAL GENERATION
    # =============================================================================
    cir_eval = cir.subs(eval_dict)
    # signal = cir_eval.kill_except('Ibulk')['S_B'].V(s)
    
    # freq_array = np.fft.fftfreq(int(1e4), 1e-4)
    
    signal_dft_dict = dict()
    for source in signal_sources:
        cir_aux = cir_eval.kill_except(source)
        signal_expr = cir_aux['S_B'].V(j*omega)(f).evalf()
        signal_lambda = sy.lambdify(f, signal_expr)
        signal_dft_dict[source] = signal_lambda(freq_array)
    
    #%%
    # =============================================================================
    # NOISE PROPAGATION
    # =============================================================================
    cir_noise = cir_eval.kill_except(*noise_sources).subs(noise_eval_dict)
    
    noise_expr_dict = dict()
    for source in noise_sources:
        cir_source = cir_noise.kill_except(source)
    
        voutput_source = cir_source['S_B'].V.n
        
        noise_expr_dict[source] = voutput_source
    
    
    # fig, ax = plt.subplots()
    
    NOISE_squared = 0 * freq_array
    for source, noise_expr in noise_expr_dict.items():
    
        # hack
        noise_fun = noise_expr(omega)(f) 
        noise_array = noise_fun.evaluate(freq_array)
        
        NOISE_squared += noise_array**2
        # ax.loglog(freq_array, noise_array, label=source)
    
    NOISE = NOISE_squared**0.5
    # ax.loglog(freq_array, NOISE, label='Total', lw=2, color='k')
    
    # ax.set_ylabel(r'LPSD [V/$\sqrt{Hz}$]')
    # ax.set_xlabel(r'Frequency [Hz]')
    # ax.legend()
    # ax.grid()
    
    #%%
    # =============================================================================
    # NEP and RESOLUTION
    # =============================================================================
    nep_dict = {k:NOISE/np.abs(v) for k,v in signal_dft_dict.items()}
    res_dict = {k:(np.trapz(4/nep**2, freq_array))**-0.5 for k,nep in nep_dict.items()}

    res_alpha.append(res_dict['Ibulk'])
