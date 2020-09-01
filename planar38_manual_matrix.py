#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 11:01:45 2020

@author: misiak
"""


import numpy as np

file_path = 'planar38_maxwell.csv'

A = np.loadtxt(file_path, comments='%', delimiter=',')

def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)

def maxwell_to_mutual(mat):
    
    ret = np.zeros(mat.shape)
    
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if i == j:
                ret[i,j] = sum(mat[i])
            else:
                ret[i,j] = -mat[i,j]
                
    return ret
            
M = maxwell_to_mutual(A)

H = np.zeros(A.shape) + 50e-12 * np.eye(A.shape[0])

AA = A + H

B = np.linalg.inv(A)

BB = np.linalg.inv(AA)

# C = B[:, 1] - B[:, 3]

# CC = BB[:, 1] - BB[:, 3]

# D = C[0] / C[1]

# DD = CC[0] / CC[1]

# print(D)
# print(DD)