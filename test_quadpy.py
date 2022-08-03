import quadpy
import sympy
from mpmath import mp
import numpy as np
# 
mp.dps = 15
scheme = quadpy.c1.gauss_legendre(3, mode="mpmath")
print(scheme.weights_symbolic)
