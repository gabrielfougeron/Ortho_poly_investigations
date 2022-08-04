import os
import quadpy
import mpmath
from mpmath import mp
import numpy as np

import orthopy

from helper import *

mp.dps = 50

max_deg = 10

# scheme_builder_list = [quadpy.c1.gauss_legendre]

scheme_builder_list = [
    quadpy.c1.gauss_legendre,
    # quadpy.c1.chebyshev_gauss_1,
    # quadpy.c1.chebyshev_gauss_2,
    ]


__PROJECT_ROOT__ = os.path.abspath(os.path.dirname(__file__))
save_dir = os.path.join(__PROJECT_ROOT__,"Butcher_tables")
if not(os.path.isdir(save_dir)):
    os.makedirs(save_dir)

for scheme_builder in scheme_builder_list:

    # for deg in [5]:
    for deg in range(1,max_deg+1):

        # integration_scheme = quadpy.c1.gauss_legendre(deg // 2 + 1,mode="mpmath")
        integration_scheme = quadpy.c1.gauss_legendre(20,mode="mpmath")
        integration_points = integration_scheme.points_symbolic
        integration_weights = integration_scheme.weights_symbolic

        scheme = scheme_builder(deg,mode="mpmath")
        print(scheme.name,deg)

        c_table = (scheme.points_symbolic + mp.mpf('1')) / mp.mpf('2')
        # b_table = scheme.weights_symbolic  / mp.mpf('2')

        bar = ComputeLagrangeBarCoeffs(c_table)
        # Lagrange_pols = lambda x : ModifiedLagrange(x,c_table,bar)
        Lagrange_pols = lambda x : BarycentricLagrange(x,c_table,bar)


        a_table = np.zeros((deg,deg),dtype=object)
        for i in range(deg):
            a_table[i,:] = integrate(Lagrange_pols,integration_weights,integration_points,a=mp.mpf('0'),b=c_table[i])

        b_table = integrate(Lagrange_pols,integration_weights,integration_points,a=mp.mpf('0'),b=mp.mpf('1'))


        # print(a_table)
        # print(b_table[0])
        # print(b_table.astype(np.float64)[0])
        # print(c_table)

        # print(np.sum(b_table))

        save_dict = {
            "a_table"       : a_table.astype(np.float64),
            "b_table"       : b_table.astype(np.float64),
            "c_table"       : c_table.astype(np.float64),
            'bar_table'     : bar.astype(np.float64)    ,
        }

        filename = os.path.join(save_dir,scheme.name.replace(" ","_")+"_"+str(deg)+'.npz')

        np.savez(filename,**save_dict)
        # np.savez_compressed(filename,**save_dict)

        print("")
