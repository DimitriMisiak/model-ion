#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 00:55:35 2020

@author: misiak
"""


import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

def custom_autoscale(axis, xdata, ydata):
    
    gold2 = ((1+5**0.5)/2)**0.5
    
    xmin, xmax = np.min(xdata), np.max(xdata)
    ymin, ymax = np.min(ydata), np.max(ydata)
    
    xcore = (xmax + xmin)/2
    ycore = (ymax + ymin)/2
    
    xside = (xmax - xmin) * gold2
    yside = (ymax - ymin) * gold2
    
    xinf, xsup = xcore - xside/2, xcore + xside/2
    yinf, ysup = ycore - yside/2, ycore + yside/2
    
    axis.set_xlim(xinf, xsup)
    axis.set_ylim(yinf, ysup)


def basic_corner(samples, labels, axes=None, num='Basic corner plot', **kwargs):
    
    assert samples.ndim == 2
    
    clen, cnum = samples.shape
    assert clen > cnum
    
    nplot = cnum - 1
    
    new_flag = False
    if axes is None:
        fig, axes = plt.subplots(nrows=nplot, ncols=nplot, figsize=(9, 9),
                                 num=num, sharex='col', sharey='row')
        new_flag = True
    
    else:
        fig = axes.flatten()[0].get_figure()
#    # lower triangle axes without diagonal
#    axes_active = np.tril(axes, -1)
#    # upper triangle axes with diagonal
#    axes_discarded = np.triu(axes, 0)
    
    options = {'ls':'none',
               'marker':'.',
               'zorder':9,
               'color':'k',
               'markersize':6,
               'alpha':0.5,
               }
    options.update(kwargs)
    
    for i in range(nplot):
        for j in range(nplot):
            ax = axes[i,j]
            
            # removing upper plots
            if i < j:
                try:
                    fig.delaxes(ax)
                except:
                    pass
        
            x_data = samples[:, (cnum-1+j)%cnum]
            y_data = samples[:, i]
            
            ax.plot(
                    x_data, y_data,
                    **options
            )            
            
            custom_autoscale(ax, x_data, y_data)
        
            ax.grid(alpha=0.3)
        
            if (i==0) and (j==0):
                ax.legend(loc='lower left', framealpha=1,
                          bbox_to_anchor=(1.05, 0.05), borderaxespad=0.,
                )
        
            if new_flag:
                if (i==nplot-1):
                    ax.set_xlabel(
                            '{}'.format(
                                    labels[(cnum-1+j)%cnum].replace('_', ' ')
                            )
                    )
                        
                if (j==0):
                    ax.set_ylabel(
                            '{}'.format(
                                    labels[i].replace('_', ' ')
                            )
                    )
    
    if new_flag:
        pass
        # fig.text(0.65, 0.98, num,
        #          horizontalalignment='center',
        #          verticalalignment='center',
        #          bbox=dict(facecolor='lime', alpha=0.5))
    
    fig.tight_layout()
    fig.subplots_adjust(hspace=.0, wspace=.0)
    
    return fig, axes

nsamples = 100

# fid38
loc_dict = {
    "bulk": np.array([-1.23, -3.52, 1.23, 3.52]),
    "veto top": np.array([1.65, -2.24, 0.11, 0.05]),
    "veto bottom": np.array([0.02, 0.01, -1.65, 2.24]),
    "equator": np.array([2.77, 1.23, -2.77, -1.23]),
}

event_dict = loc_dict.copy()

# # fid803
# loc_dict = {
#     "bulk": np.array([-0.47337116, -0.73655778,  0.47358752,  0.73505918]),
#     "veto top": np.array([ 0.17285785, -0.25659646,  0.00833141,  0.00659015]),
#     "veto bottom": np.array([-0.00759521, -0.00585395, -0.17343654,  0.25561771]),
#     "equator": np.array([ 0.6386338 ,  0.47410736, -0.63869265, -0.47285132]),
# }

# mu_dict = {k:abs(0.1*v) for k,v in loc_dict.items()}

samples_dict = dict()
for k in loc_dict.keys():
    loc_array = loc_dict[k]
    # mu_array = mu_dict[k]
    mu_array = 0.001
    samples_dict[k] = np.random.normal(loc=loc_array, scale=mu_array, size=(nsamples, 4))

fig, axes = basic_corner(samples_dict['bulk'], labels='ABCD', color='limegreen', marker='s', label='Bulkw/ 0pF')
fig, axes = basic_corner(samples_dict['veto top'], axes=axes, labels='ABCD', color='orange',marker='s', label='Veto Topw/ 0pF')
fig, axes = basic_corner(samples_dict['veto bottom'], axes=axes, labels='ABCD', color='crimson',marker='s', label='Veto Bottomw/ 0pF')
fig, axes = basic_corner(samples_dict['equator'], axes=axes, labels='ABCD', color='slateblue',marker='s', label='Equatorw/ 0pF')

for ax in np.ravel(axes):
    ax.axvline(0, ls='--', color='k')
    ax.axhline(0, ls='--', color='k')

# fid38
loc_dict = {
    "bulk": np.array([-0.05, -0.80, 0.05, 0.80]),
    "veto top": np.array([0.69, -0.75, 0, 0]),
    "veto bottom": np.array([0, 0, -0.69, 0.75]),
    "equator": np.array([0.75, 0.05, -0.75, -0.05]),
}
    
# # fid803
# loc_dict = {
#     "bulk": np.array([-0.06945887, -0.20961447,  0.06958691,  0.20946544]),
#     "veto top": np.array([ 0.12903953, -0.13797299,  0.00336613,  0.00218261]),
#     "veto bottom": np.array([-0.00325698, -0.00207347, -0.12908928,  0.13780506]),
#     "equator": np.array([ 0.19524142,  0.06956801, -0.19531007, -0.06947776]),
# }
# mu_dict = {k:abs(0.1*v) for k,v in loc_dict.items()}

color_dict = {
    "bulk": 'limegreen',
    "veto top": 'orange',
    "veto bottom": 'crimson',
    "equator": 'slateblue',
}
samples_dict = dict()
for k in loc_dict.keys():
    loc_array = loc_dict[k]
    # mu_array = mu_dict[k]
    mu_array = 0.001
    samples_dict[k] = np.random.normal(loc=loc_array, scale=mu_array, size=(nsamples, 4))

fig, axes = basic_corner(samples_dict['bulk'], axes=axes, labels='ABCD', color='limegreen', marker='o',mec='k',  label='Bulk w/ 150pF')
fig, axes = basic_corner(samples_dict['veto top'], axes=axes, labels='ABCD', color='orange',marker='o', mec='k', label='Veto Top w/ 150pF')
fig, axes = basic_corner(samples_dict['veto bottom'], axes=axes, labels='ABCD', color='crimson',marker='o', mec='k', label='Veto Bottom w/ 150pF')
fig, axes = basic_corner(samples_dict['equator'], axes=axes, labels='ABCD', color='slateblue',marker='o', mec='k', label='Equator w/ 150pF')


# fig, axes = basic_corner(samples_dict['bulk'], axes=axes, labels='ABCD', color='k', label='w/ 50pF')
# fig, axes = basic_corner(samples_dict['veto top'], axes=axes, labels='ABCD', color='k')
# fig, axes = basic_corner(samples_dict['veto bottom'], axes=axes, labels='ABCD', color='k')
# fig, axes = basic_corner(samples_dict['equator'], axes=axes, labels='ABCD', color='k')



### Bonus

maxwell_matrix = np.loadtxt('fid38_maxwell.txt', comments='%')

A = np.dot(maxwell_matrix, event_dict['bulk'])