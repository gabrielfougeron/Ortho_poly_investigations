from mpmath import mp

# weights and points are assumed to be defined on (-1,1)
def integrate(fun,weights,points,a=mp.mpf('-1'),b=mp.mpf('1')):
    
    mid = (a+b) / mp.mpf('2')
    stretch = (b-a) / mp.mpf('2')

    new_points = [ mid + stretch * x for x in points]
    new_weights = [ stretch * w for w  in weights]

    fun_eval = [fun(x) for x in new_points]

    return sum( [w*f for (w, f) in zip(new_weights, fun_eval)])


def ComputeLagrangeBarCoeffs(points):

    bar = []
    n = len(points)

    for i in range(len(n)):

        coeff = mp.mpf('1')

        for j in range(len(n)):

            if (i != j):
                coeff = coeff * (points[i] - points(j))

        bar.append( mp.mpf('1') / coeff)

    return bar

def ModifiedLagrange(x,points,bar=None):
    if bar is None:
        bar = ComputeLagrangeBarCoeffs(points)
    
    n = len(points)

    l = mp.mpf('1')
    p = []
    
    for i in range(len(n)):

        dx = (x - point[i])
        l = l * dx
        p.append(bar/dx)    

    for i in range(len(n)):
        
        p[i] = p[i] * l

    return p

def BarycentricLagrange(x,points,bar=None):
    if bar is None:
        bar = ComputeLagrangeBarCoeffs(points)

    n = len(points)

    l = mp.mpf('1')
    p = []
    
    for i in range(len(n)):

        p.append(bar/dx)   

    the_sum = mp.mpf('0')

    for i in range(len(n)):

        the_sum = the_sum + p[i]

    for i in range(len(n)):
        
        p[i] = p[i] / the_sum

    return p

