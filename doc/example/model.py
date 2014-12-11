def density(x,y,z):
    import math
    r = math.sqrt(x*x+y*y+z*z)
    c = 1.5e6*math.pow(r/(300*1.49598e11),-1.5)*1e6
    return c

def temperature(x,y,z):
    import numpy
    import math
    x0 = 0
    temp=numpy.array([[2.0e13, 5.0e13, 8.0e13, 1.1e14, 1.4e14, 1.7e14, 2.0e14, 2.3e14, 2.6e14, 2.9e14],[44.777, 31.037, 25.718, 22.642, 20.560, 19.023, 17.826, 16.857, 16.050, 15.364]])
    r = math.sqrt(x*x+y*y+z*z)
    if ( r > temp[0][0]) and ( r < temp[0][9] ):
        for i in range(0, 9):
            if (r > temp[0][i] ) and (r < temp[0][i+1] ):
                x0=i
                break
    if r<temp[0][0]:
        return temp[1][0]
    elif r>temp[0][9]:
        return temp[1][9]
    else:
        return (temp[1][x0]+(r-temp[0][x0])*(temp[1][x0+1]-temp[1][x0])/(temp[0][x0+1]-temp[0][x0]))

def velocity(x,y,z):
    import math
    R = math.sqrt(x*x+y*y+z*z)
    theta = math.atan2(math.sqrt(x*x+y*y), z)
    phi=math.atan2(y,x)
    r=-math.sqrt(2*6.67e-11*1.989e30/R)
    return ( r*math.sin(theta)*math.cos(phi), r*math.sin(theta)*math.sin(phi), r*math.cos(theta) )

def doppler(x,y,z):
    return 200

def abundance(x,y,z):
    return 1.e-9
