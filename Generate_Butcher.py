import quadpy
from mpmath import mp
import numpy as np

import orthopy

from helper import *


mp.dps = 100

max_deg = 5

scheme_builder_list = [quadpy.c1.gauss_legendre]

for scheme_builder in scheme_builder_list:

    for deg in range(1,max_deg):

        integration_scheme = quadpy.c1.gauss_legendre(deg // 2 + 1,mode="mpmath")
        integration_points = integration_scheme.points_symbolic
        integration_weights = integration_scheme.weights_symbolic

        scheme = scheme_builder(deg,mode="mpmath")
        print(deg)
        print(scheme.name)

        c_table = (scheme.points_symbolic + mp.mpf('1')) / mp.mpf('2')
        b_table = scheme.weights_symbolic  / mp.mpf('2')

        bar = ComputeLagrangeBarCoeffs(c_table)
        Lagrange_pols = lambda x : ModifiedLagrange(x,c_table,bar)
        # Lagrange_pols = lambda x : BarycentricLagrange(x,c_table,bar)


        a_table = np.zeros((deg,deg),dtype=object)
        for i in range(deg):
            a_table[i,:] = integrate(Lagrange_pols,integration_weights,integration_points,a=mp.mpf('0'),b=c_table[i])


        # print(a_table)
        # print(b_table)
        # print(c_table)

        # print(type(a_table))
        # print(type(b_table))
        # print(type(c_table))


        print("")
