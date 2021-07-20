#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:38:58 2020

@author: misiak
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


plt.close('all')

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1

# # =============================================================================
# # Extra detector capacitance
# # =============================================================================

# capa_array = np.array([3, 8, 13, 28])
# res_array = np.array([34.6, 32.9, 34, 42])

# plt.figure()
# plt.plot(capa_array, res_array, label='B', marker='o')

# plt.ylabel('Resolution [eV]')
# plt.xlabel('Hemt Line Equivalent Capacitance [pF]')
# plt.grid()


# # =============================================================================
# # ALPHA multiplicative factor
# # =============================================================================
# alpha_array = np.array([0.25, 0.5, 0.75, 1, 1.5, 2])
# res_array = np.array([19, 24, 29, 34, 44, 53])

# plt.figure()
# plt.plot(alpha_array, res_array, label='B', marker='o')

# plt.ylabel('Resolution [eV]')
# plt.xlabel('Capacitance matrix Multiplier Factor')

# plt.ylabel('Resolution [eV]')
# plt.xlabel('Multiplicative factor')
# plt.grid()



# # =============================================================================
# # BETA multiplicative factor
# # =============================================================================
# alpha_array = np.array([0.25, 0.5, 0.75, 1, 1.5, 2])
# res_array = np.array([20.6, 25, 29.5, 34, 43, 52.2])

# plt.figure()
# plt.plot(alpha_array, res_array, label='B', marker='o')

# plt.ylabel('Resolution [eV]')
# plt.xlabel('Capacitance matrix Multiplier Factor')

# plt.ylabel('Resolution [eV]')
# plt.xlabel('Multiplicative factor')
# plt.grid()

# =============================================================================
# STACKED PLOTS
# =============================================================================

### FIGURE
fig, axes = plt.subplots(
    figsize=(6.3,3.9),
    ncols=2,
    sharey=True
)

ax0, ax1 = axes

capa_array = np.array([0, 2.5, 5, 7.5, 10, 12.25, 15, 17.25, 20, 25, 30])
res_array = 1000 * np.array([0.028687532243911417, 0.028666816576754683, 0.029120845465633732, 0.02985781984155704, 0.030775651830789987, 0.031707279733052485, 0.03294183393358135, 0.03400988762469854, 0.0353681711132052, 0.03794512441889855, 0.04061757019586626])



# ([0.03290383600972107, 0.033202360410120826, 0.03397135433644549, 0.035016803102193524, 0.03623727033358495, 0.03743724039733141, 0.03899577887540935, 0.040326420506200555, 0.042003390354634886, 0.04515548082984242, 0.048399821174430646])
ax0.plot(capa_array, res_array, label='B', marker='.', color='slateblue')
ax0.plot(capa_array, res_array/2**0.5, label='B', marker='.', color='coral')
# capa_array = np.array([0, 2.5, 5, 7.5, 10, 12.25, 15])
# res_array = 1000 * np.array([0.03287831840420088, 0.033188303577250215, 0.03397135433644549, 0.035032416003738236, 0.03626949566219288, 0.037485004836537716, 0.03906307631620481])
# ax0.plot(capa_array, res_array, label='B', marker='.')

alpha_array = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2])
# res_array = np.array([
#     0.016018174543648886,
#     0.01802749073727444,
#     0.020043048368436348,
#     0.02205711920992607,
#     0.024065434874996485,
#     0.02606568382870057,
#     0.02805669717907899,
#     0.030037988252031807,
#     0.032009482520391784,
#     0.03397135433644549,
#     0.03883624633064418,
#     0.04364983860477241,
#     0.0484193178583591,
#     0.05315133586812252
# ]) * 1000

res_array = np.array([
    0.014097173621937239,
    0.015778640734432494,
    0.01746681746401001,
    0.01915412507044838,
    0.020836324369383444,
    0.022511093710891526,
    0.02417724444713354,
    0.025834270938654265,
    0.02748208351811249,
    0.029120845465633732,
    0.0331804196923417,
    0.03719215635608707,
    0.04116320033635307,
    0.045100152360066564
]) * 1000

ax1.plot(alpha_array, res_array, label='B', marker='.', color='slateblue')
ax1.plot(alpha_array, res_array/2**0.5, label='B', marker='.', color='coral')

### Axes labels
ax0.set_ylabel(r'Ionization Energy Resolution $\sigma$ / eV')
ax0.set_xlabel(r'Cabling Capacitance / pF')
ax1.set_xlabel(r'Multiplicative factor $\alpha$')

ax0.spines['right'].set_linewidth(2)

ax0.yaxis.set_major_locator(mticker.MultipleLocator(10))
ax0.yaxis.set_minor_locator(mticker.MultipleLocator(2.5))

ax0.xaxis.set_major_locator(mticker.MultipleLocator(10))
ax0.xaxis.set_minor_locator(mticker.MultipleLocator(2.5))

ax1.xaxis.set_major_locator(mticker.MultipleLocator(0.5))
ax1.xaxis.set_minor_locator(mticker.MultipleLocator(.25))

# ### Axes scale
# ax0.set_yscale('log')
# ax0.set_xscale('log')

# ax0.set_xlim(freq_array.min(), freq_array.max())
# ax0.set_ylim(1e-11, 3e-7)
# ax_twin.set_ylim(1e-11, 3e-7)

# # subtitles in legend
# import matplotlib.text as mtext
# class LegendTitle(object):
#     def __init__(self, text_props=None):
#         self.text_props = text_props or {}
#         super(LegendTitle, self).__init__()

#     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#         x0, y0 = handlebox.xdescent, handlebox.ydescent
#         # title = mtext.Text(x0, y0, r'\underline{' + orig_handle + '}',
#         #                    usetex=True, **self.text_props)
#         title = mtext.Text(x0, y0, orig_handle,
#                            usetex=True, **self.text_props)
#         handlebox.add_artist(title)
#         return title


### Legends
ax0.legend(
    handles=[
        plt.Line2D([], [], lw=1, color='slateblue', marker='.', ls='-'),
        plt.Line2D([], [], lw=1, color='coral', ls='-', marker='.'),
    ],
    labels=[
        'Single Collect Channel Resolution $\sigma (B)$',
        r'Combined Collect Channel Resolution $\sigma (\frac{B+D}{2})$'
    ],
    loc='lower center', bbox_to_anchor=(1, 1), frameon=False
)

### Grid
for ax in fig.axes:
    ax.grid(True, alpha=0.5, which='major')
    ax.grid(True, alpha=0.1, which='minor')

### Figure adjustments
fig.align_ylabels(fig.axes)    
fig.tight_layout()
fig.subplots_adjust(wspace=.0)

### Saving
fig.savefig('pl38_scanplot.pdf')




