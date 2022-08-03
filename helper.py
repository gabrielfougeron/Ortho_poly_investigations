from mpmath import mp
import numpy as np

# weights and points are assumed to be defined on (-1,1)
def integrate(fun,weights,points,a=mp.mpf('-1'),b=mp.mpf('1')):
    
    n = len(points)

    mid = (a+b) / mp.mpf('2')
    stretch = (b-a) / mp.mpf('2')

    for i in range(n):

        x = mid + stretch * points[i]
        w = stretch * weights[i]

        f = fun(x)

        if (i==0):
                
            the_quad = np.array([ mp.mpf('0') for j in range(len(f))],dtype=object)

        for j in range(len(f)):

            the_quad[j] = the_quad[j] + w * f[j]

    return the_quad


def ComputeLagrangeBarCoeffs(points):

    n = len(points)
    bar = np.array([ mp.mpf('0') for j in range(n)],dtype=object)

    for i in range(n):

        coeff = mp.mpf('1')

        for j in range(n):

            if (i != j):
                coeff = coeff * (points[i] - points[j])

        bar[i] =  mp.mpf('1') / coeff

    return bar

def ModifiedLagrange(x,points,bar=None):
    if bar is None:
        bar = ComputeLagrangeBarCoeffs(points)
    
    n = len(points)

    l = mp.mpf('1')
    p = np.array([ mp.mpf('0') for j in range(n)],dtype=object)

    try:
        
        for i in range(n):

            dx = (x - points[i])
            l = l * dx
            p[i] = bar[i] / dx

        for i in range(n):
            
            p[i] = p[i] * l

    except ZeroDivisionError: # occured at index i

        p[i]= mp.mpf('1')

    return p



def BarycentricLagrange(x,points,bar=None):
    if bar is None:
        bar = ComputeLagrangeBarCoeffs(points)

    n = len(points)

    l = mp.mpf('1')
    p = np.array([ mp.mpf('0') for j in range(n)],dtype=object)
    
    try:
            
        for i in range(n):
            dx = (x - points[i])
            p[i] = bar[i] / dx

        the_sum = mp.mpf('0')

        for i in range(n):

            the_sum = the_sum + p[i]

        for i in range(n):
            
            p[i] = p[i] / the_sum

    except ZeroDivisionError: # occured at index i

        p[i]= mp.mpf('1')


    return p
