#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:56:53 2020

@author: misiak

For the 3d Wheastone bridge with capacitor:
    
see notebook for illustration
"""


import sympy as sy

# impedance
z_list = sy.symbols(
    'z11, z12, z13, z14, z22, z23, z24, z33, z34, z44'
)
z11, z12, z13, z14, z22, z23, z24, z33, z34, z44 = z_list

# admittance
y_list = sy.symbols(
    'y11, y12, y13, y14, y22, y23, y24, y33, y34, y44'
)
y11, y12, y13, y14, y22, y23, y24, y33, y34, y44 = y_list

# capacitance
c_list = sy.symbols(
    'c11, c12, c13, c14, c22, c23, c24, c33, c34, c44'
)
c11, c12, c13, c14, c22, c23, c24, c33, c34, c44 = c_list

w = sy.symbols('w')
subs_impedance = {z: 1/(sy.I * w * c) for z,c in zip(z_list, c_list)}
subs_admittance = {y: sy.I * w * c for y,c in zip(y_list, c_list)}

v0, v1, v2, v3, v4 = sy.symbols('v:5')

subs_dict = {
    'v1': 1, #V
    'v0': 0, #V
}


# using impedance
subs_aux = subs_impedance
# Kirchhoff's current law in fourier space
# at node 2:
eq_node2 = (v1-v2)/z12 + (v3-v2)/z23 + (v4-v2)/z24 + (v0-v2)/z22
# # at note 3:
eq_node3 = (v1-v3)/z13 + (v2-v3)/z23 + (v4-v3)/z34 + (v0-v3)/z33
# # at note 3:
eq_node4 = (v1-v4)/z14 + (v2-v4)/z24 + (v3-v4)/z34 + (v0-v4)/z44
# total current
Y = (v2-v0)/z22 + (v3-v0)/z33 + (v4-v0)/z44

# # # using admittance instead of impedance, for alternative expression
# subs_aux = subs_admittance
# # at node 2:
# eq_node2 = (v1-v2)*y12 + (v3-v2)*y23 + (v4-v2)*y24 + (v0-v2)*y22
# # # at note 3:
# eq_node3 = (v1-v3)*y13 + (v2-v3)*y23 + (v4-v3)*y34 + (v0-v3)*y33
# # # at note 3:
# eq_node4 = (v1-v4)*y14 + (v2-v4)*y24 + (v3-v4)*y34 + (v0-v4)*y44
# # total current
# Y = (v2-v0)*y22 + (v3-v0)*y33 + (v4-v0)*y44

# solving the system
solution_dict = sy.solve([eq_node2, eq_node3, eq_node4], [v2, v3, v4])

#equivalent resistance
Y_eval = (Y.subs(solution_dict)).subs(subs_dict).simplify()
Z_eval = 1/Y_eval

Y_aux = Y_eval.subs(subs_aux).simplify()
Z_aux = 1/Y_aux

# print('R=')
# sy.pprint(R_eval.simplify())