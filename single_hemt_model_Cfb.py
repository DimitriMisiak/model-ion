#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:59:15 2020

@author: filippini

Simple Hemt model with a single electrode.
The Hemt has a resistive feedback.
"""

import lcapy
from lcapy import Circuit,f, j, t, omega
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

freq_array = 10**np.linspace(0, 4, int(1e3))

eval_dict = {
    'Cd': 14.926123e-12, #F 
    'Cp': 5e-12, #F
    'Rpolar': 1e10, #Ohms
    'Cc': 2e-9, #F
    'Chemt': 5e-12, #F
    'Rfb': 1e10, #1e10, #F
    'Cfb': 3e-12, #F
    'qe': 1.60e-19, #C,
    'Npair': 1e3/3, # number of pairs in 1keV
    'dt': 1e-5
}

cir = Circuit("""
    O S 0; down
    W 0 0_g; down=0.2, cground                   

    W E 1_d; right
    Cd 1_d 0_d; down
    W 1_d 1_p; right
    Cp 1_p 0_p; down
    W 1_p 1_r; right
    Rpolar 1_r 2_r; down
    Vj 2_r 3_r noise; down, l=e_J, color=red
    Vdaq 3_r 0_r noise; down, l=e_{DAQ}, color=red
    
    Cc 1_r 1_fb; right=2
    Cfb 1_fb 0_fb; down
    
    W 1_fb 1_h; right
    Chemt 1_h 0_h; down
    W 1_h 1_s; right

    Vn 1_s S noise; right, color=red, l=e_n
    In 1_s 0_s noise; down, color=red, l=i_n
    
    W 0_E 0_d; right
    W 0_d 0_p; right
    W 0_p 0_r; right
    W 0_r 0_fb; right
    W 0_fb 0_h; right
    W 0_h 0_s; right
    W 0_s 0; right
    
    Iq E 0_E {Npair*qe*delta(t)}; down, color=blue, l_=I_{q}
    
    ; draw_nodes=primary
    ; label_nodes=alpha
    ; style=european
    ; label_ids=False
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

# # Feedback impedance Zfb is a resistance
# cir_rfb = cir.subs({'Zfb':'Rfb'})

# # Feedback impedance Zfb is a capacitive relay
# cir_cfb = cir.subs({'Zfb':'(j*omega*Cfb)**-1'}).kill('Vfb')

# circuit without source
cir_kill = cir.kill()

#%%
# =============================================================================
# SIGNAL GENERATION
# =============================================================================
cir_signal = cir.kill_except('Iq').subs(eval_dict)

# temporal plot
time_array = np.linspace(-0.5, 1.5, int(1e4))
iinput_array_temp = cir_signal.Iq.I(t).evaluate(time_array)
vinput_array_temp = cir_signal['E'].V(t).evaluate(time_array)
voutput_array_temp = cir_signal['S'].V(t).evaluate(time_array)

fig, ax = plt.subplots(nrows=2, sharex=True)
ax[0].plot(time_array, iinput_array_temp, label='Current Input from charge collection')
ax[0].set_ylabel('Current [A]')
ax[1].plot(time_array, vinput_array_temp, label='Input, Node E')
ax[1].plot(time_array, voutput_array_temp, label='Output, Node S')
ax[1].set_ylabel('Voltage [V]')
ax[1].set_xlabel('Time [s]')

for a in ax:
    a.grid()
    a.legend()

# frequency plot
iinput_array_freq = abs(cir_signal.Iq.I(j*omega)(f).evaluate(freq_array))
vinput_array_freq = abs(cir_signal['E'].V(j*omega)(f).evaluate(freq_array))
voutput_array_freq = abs(cir_signal['S'].V(j*omega)(f).evaluate(freq_array))

SIGNAL = voutput_array_freq

fig, ax = plt.subplots(nrows=2, sharex=True)
ax[0].plot(freq_array, iinput_array_freq, label='Current Input from charge collection')
ax[0].set_ylabel('Current [A]')
ax[1].loglog(freq_array, vinput_array_freq, label='Input, Node E')
ax[1].loglog(freq_array, voutput_array_freq, label='Output, Node S')
ax[1].set_ylabel('Voltage [V]')
ax[1].set_xlabel('Frequency [Hz]')

for a in ax:
    a.grid()
    a.legend()

#%%
# =============================================================================
# NOISE PROPAGATION
# =============================================================================

cir_noise_symbol = cir.kill('Iq')
cir_noise = cir_noise_symbol.subs(eval_dict)

tot_lab = 'Tot'
noise_expr_dict = dict()
for source in (cir_noise.sources + [tot_lab,]) :
    
    if source == tot_lab:
        cir_source = cir_noise
    else:
        cir_source = cir_noise.kill_except(source)

    voutput_source = cir_source['S'].V.n
    
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


#%%
# =============================================================================
# NEP and RESOLUTION
# =============================================================================
NEP = NOISE / SIGNAL
invnep2 = 4/NEP**2
RES = (np.trapz(invnep2, freq_array))**-0.5

print('Resolution = ', RES*1e3, ' eV')

### PLOT
noise_color = 'coral'
signal_color = 'slateblue'

fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(10,10))

ax = axes[0]
ax.tick_params(axis='y', labelcolor=noise_color)
ax.loglog(freq_array, NOISE, color=noise_color, lw=2)
ax.set_ylabel('Noise LPSD [V/$\sqrt{Hz}$]', color=noise_color)
ax_twin = ax.twinx()
ax_twin.tick_params(axis='y', labelcolor=signal_color)
ax_twin.loglog(freq_array, SIGNAL, color=signal_color, lw=2)
ax_twin.set_ylabel('Signal [V]', color=signal_color)

# ax.set_xlabel('Frequency [Hz]')
ax.grid()

ax = axes[1]
ax.loglog(freq_array, NEP, color='crimson')
ax.set_ylabel('Noise Equivalent Power (NEP)')
ax.set_xlabel('Frequency [Hz]')
ax.grid()

axes[0].set_title(
    'Resolution Calculation\n'
    + r'$\sigma_{{ion}}$ = {} eV'.format(RES*1e3)
)

fig.tight_layout()
fig.subplots_adjust(hspace=.0)
