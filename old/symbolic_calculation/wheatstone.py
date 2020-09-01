#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:56:53 2020

@author: misiak

For the Wheastone bridge:
    
           2
    +-+R1+---+R4+-+
    |      |      |
    |      +      |
1+--+      R3     +--+4
    |      +      |
    |      |      |
    +-+R2+---+R5+-+
           3

This calculation is correct :)
"""


import sympy as sy

r1, r2, r3, r4, r5 = sy.symbols('r1:6')
v1, v2, v3, v4 = sy.symbols('v1:5')



subs_dict = {
    'v1': 1, #V
    'v4': 0, #V
}

# Kirchhoff's current law
# at node 2:
eq_node2 = (v1-v2)/r1 + (v4-v2)/r4 + (v3-v2)/r3
# at note 3:
eq_node3 = (v1-v3)/r2 + (v4-v3)/r5 + (v2-v3)/r3

# using sigma instead of r, for easy crosschek
# eq_node2 = (v1-v2)*r1 + (v4-v2)*r4 + (v3-v2)*r3
# eq_node3 = (v1-v3)*r2 + (v4-v3)*r5 + (v2-v3)*r3


# total current
I = (v3-v4)*r5 + (v2-v4)*r4
#equivalent resistance
R = 1/I

solution_dict = sy.solve([eq_node2, eq_node3], [v2, v3])

R_eval = (R.subs(solution_dict)).subs(subs_dict)

print('R=')
sy.pprint(R_eval.simplify())