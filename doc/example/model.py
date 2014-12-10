def density(x,y,z):
    import math
    r = math.sqrt(x*x+y*y+z*z)
    c = 1.5e6*math.pow(r/(300*1.49598e11),-1.5)*1e6
    return c

def temperature(a,b):
    c = a*b
    return c

def velocity(a,b):
    c = a*b
    return c

def doppler(a,b):
    c = a*b
    return c

def abundance(a,b):
    c = a*b
    return c
