#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:59:15 2020

@author: filippini
"""

import lcapy
import numpy as np
import matplotlib.pyplot as plt
from ionisation_elcapy_fct import resolution



#noise hemt 
En = lcapy.expr('(e0**2 +(ea/f**0.5)**2 + (eb/f)**2)**0.5 ')
In = lcapy.expr('(i0**2 +(ia*(f)**0.5)**2 + (ib*f)**2)**0.5 ')
e_dict = {'e0': 0.22e-9, 'ea': 36e-9, 'eb': 0}
i_dict = {'i0': 4.0e-5*1e-18, 'ia':2.6e-18, 'ib':0}
En_eval = En.subs(e_dict)
In_eval = In.subs(i_dict)
en_noise = En_eval(lcapy.f)
in_noise = In_eval(lcapy.f)



Cd = 10e-12
Cp = 10e-12
Rpolar = 1e10
Cc = 2e-9
Chemt = 5e-12
Cfb =  1e-17
T = 20e-3
k = 1.38e-23
e = 1.6e-19
tension = 333*e


frequency = np.arange(1,50e3,1)


circuit = lcapy.Circuit("""
                        Cd 1 0; down
                        W 1 1_2; right
                        Cp 1_2 0_1; down
                        W 1_2 1_3; right
                        Rpolar 1_3 0_2; down
                        Cc 1_3 2; right
                        W 2 2_1; right

                        W 2_1 2_2; rigth
                        Chemt 2_2 0_4; down
                        W 2_2 2_3; right
                        Vn 2_3 3 noise; right
                        In 2_3 0_5 noise; down
                        W 0 0_1; right
                        W 0_1 0_2; right
                        W 0_2 0_3; right
                        W 0_3 0_4; right
                        W 0_4 0_5; right
                        """)


# circuit.draw()
R_noise = circuit.noisy()
R_noise.draw()

print('a')
R_noise = R_noise.subs({'Cd':Cd, 'Cp':Cp, 'Rpolar':Rpolar,
                          'Cc':Cc, 'Chemt':Chemt,
                          # 'Cfb':Cfb,
                          'T':T,
                          'k':k,
                          })

#calcul impedance
impedance = R_noise.impedance(1,0)(lcapy.f).evaluate(frequency)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(frequency, abs(impedance), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)
print('b')

#%%calcul total noise
fig = plt.figure()
ax = fig.add_subplot(111)
for sources in R_noise.sources:
    Noise = R_noise.kill_except(sources)[3].V.n
    Noise = Noise.subs({'Cd':Cd, 'Cp':Cp, 'Rpolar':Rpolar,
                        'Cc':Cc, 'Chemt':Chemt,
                        # 'Cfb':Cfb, 
                        'T':T,
                        'k':k,
                        'Vn':en_noise(lcapy.omega),
                        'In':in_noise(lcapy.omega)})


    noise = Noise(lcapy.f).evaluate(frequency)
    ax.loglog(frequency, noise, label = f'{sources}')

Noise = R_noise[3].V.n.subs({'Cd':Cd, 'Cp':Cp, 'Rpolar':Rpolar,
                             'Cc':Cc, 'Chemt':Chemt,
                             # 'Cfb':Cfb, 
                             'T':T,
                             'k':k,
                             'Vn':en_noise(lcapy.omega),
                          'In':in_noise(lcapy.omega)})
noise = Noise(lcapy.f).evaluate(frequency)
ax.loglog(frequency, noise, label = f'Total noise')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Noise (V/sqrt(Hz))')
ax.legend()

from lcapy import C,Z,R
imped = lcapy.LSection((Z(Cc)),(Z(Chemt)))
imped_f_domain = imped.Vtransfer.frequency_response()

# on calcul pour toutes les freqs
imped_f_domain = float(R_noise.transfer(1,0,2,0).evalf())*impedance
print('res', resolution(noise, imped_f_domain, frequency))
# vf = circuit.R.v.evaluate(frequency)