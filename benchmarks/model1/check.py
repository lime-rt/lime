#!/usr/bin/env python
# check.py -- Check the results of the benchmark model
#
# This file is part of LIME, the versatile line modeling engine
#
#  Copyright (C) 2006-2014 Christian Brinch
#  Copyright (C) 2015-2016 The LIME development team


import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 14})

def readpop(filename):
    """
    Return the fractional populations a a function of position
    """

    data = np.loadtxt(filename)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    radius *= 1e2 # m -> cm
    pop = data[:, 7:]

    return radius, pop

def main():

    # Plot the LIME results

    filename1a = "model1a.pop"
    radius1a, pop1a = readpop(filename1a)

    filename1b = "model1b.pop"
    radius1b, pop1b = readpop(filename1b)

    plt.figure(figsize=(10, 9))
    plt.title("Benchmark model 1")
    plt.xlabel("r (cm)")
    plt.ylabel("Fractional population of the upper level")
    plt.xlim(1e15, 1e19)
    plt.ylim(.1, 1)
    plt.yscale('log')
    plt.xscale('log')
    p1 = plt.plot(radius1a, pop1a[:, 1], '.', label="Model 1a (LIME)")
    p2 = plt.plot(radius1b, pop1b[:, 1], '.', label="Model 1b (LIME)")

    # Overplot the results from Fig 2 of van Zadelhoff (2002) paper.
    # Values are obtained with WebPlotDigitizer.

    fiducial1a = np.array([
        [1.41e+15, 6.62e-1],
        [1.68e+15, 6.62e-1],
        [2.01e+15, 6.62e-1],
        [2.40e+15, 6.62e-1],
        [2.86e+15, 6.62e-1],
        [3.42e+15, 6.61e-1],
        [4.08e+15, 6.60e-1],
        [4.88e+15, 6.58e-1],
        [5.82e+15, 6.55e-1],
        [6.95e+15, 6.49e-1],
        [8.30e+15, 6.41e-1],
        [9.92e+15, 6.30e-1],
        [1.18e+16, 6.15e-1],
        [1.41e+16, 5.95e-1],
        [1.69e+16, 5.71e-1],
        [1.99e+16, 5.42e-1],
        [2.29e+16, 5.15e-1],
        [2.58e+16, 4.89e-1],
        [2.86e+16, 4.64e-1],
        [3.16e+16, 4.39e-1],
        [3.45e+16, 4.16e-1],
        [3.74e+16, 3.94e-1],
        [4.03e+16, 3.73e-1],
        [4.33e+16, 3.53e-1],
        [4.64e+16, 3.34e-1],
        [4.96e+16, 3.16e-1],
        [5.29e+16, 2.99e-1],
        [5.63e+16, 2.83e-1],
        [6.00e+16, 2.67e-1],
        [6.37e+16, 2.54e-1],
        [6.76e+16, 2.41e-1],
        [7.17e+16, 2.29e-1],
        [7.64e+16, 2.16e-1],
        [8.14e+16, 2.05e-1],
        [8.67e+16, 1.95e-1],
        [9.32e+16, 1.84e-1],
        [1.01e+17, 1.74e-1],
        [1.09e+17, 1.65e-1],
        [1.19e+17, 1.57e-1],
        [1.31e+17, 1.48e-1],
        [1.46e+17, 1.41e-1],
        [1.67e+17, 1.34e-1],
        [1.96e+17, 1.27e-1],
        [2.34e+17, 1.22e-1],
        [2.80e+17, 1.19e-1],
        [3.34e+17, 1.17e-1],
        [3.99e+17, 1.16e-1],
        [4.76e+17, 1.15e-1],
        [5.69e+17, 1.14e-1],
        [6.79e+17, 1.13e-1],
        [8.11e+17, 1.13e-1],
        [9.69e+17, 1.13e-1],
        [1.16e+18, 1.13e-1],
        [1.38e+18, 1.13e-1],
        [1.65e+18, 1.13e-1],
        [1.97e+18, 1.13e-1],
        [2.35e+18, 1.13e-1],
        [2.81e+18, 1.13e-1],
        [3.36e+18, 1.13e-1],
        [4.01e+18, 1.13e-1],
        [4.79e+18, 1.13e-1],
        [5.72e+18, 1.13e-1],
        [6.83e+18, 1.13e-1],
        [7.62e+18, 1.13e-1]])

    fiducial1b = np.array([
        [1.37e+15, 6.61e-1],
        [1.64e+15, 6.61e-1],
        [1.96e+15, 6.61e-1],
        [2.34e+15, 6.61e-1],
        [2.79e+15, 6.61e-1],
        [3.34e+15, 6.61e-1],
        [3.99e+15, 6.61e-1],
        [4.76e+15, 6.61e-1],
        [5.68e+15, 6.61e-1],
        [6.79e+15, 6.61e-1],
        [8.11e+15, 6.61e-1],
        [9.68e+15, 6.61e-1],
        [1.16e+16, 6.61e-1],
        [1.38e+16, 6.61e-1],
        [1.65e+16, 6.61e-1],
        [1.97e+16, 6.60e-1],
        [2.35e+16, 6.59e-1],
        [2.81e+16, 6.56e-1],
        [3.36e+16, 6.52e-1],
        [4.01e+16, 6.46e-1],
        [4.79e+16, 6.36e-1],
        [5.72e+16, 6.24e-1],
        [6.83e+16, 6.07e-1],
        [8.15e+16, 5.84e-1],
        [9.70e+16, 5.58e-1],
        [1.14e+17, 5.29e-1],
        [1.31e+17, 5.00e-1],
        [1.46e+17, 4.75e-1],
        [1.62e+17, 4.52e-1],
        [1.79e+17, 4.24e-1],
        [1.95e+17, 4.03e-1],
        [2.11e+17, 3.81e-1],
        [2.29e+17, 3.58e-1],
        [2.44e+17, 3.43e-1],
        [2.59e+17, 3.24e-1],
        [2.72e+17, 3.07e-1],
        [2.93e+17, 2.91e-1],
        [3.10e+17, 2.76e-1],
        [3.32e+17, 2.59e-1],
        [3.59e+17, 2.43e-1],
        [3.80e+17, 2.31e-1],
        [4.04e+17, 2.18e-1],
        [4.30e+17, 2.06e-1],
        [4.59e+17, 1.95e-1],
        [4.93e+17, 1.85e-1],
        [5.30e+17, 1.75e-1],
        [5.72e+17, 1.66e-1],
        [6.23e+17, 1.58e-1],
        [6.97e+17, 1.48e-1],
        [8.06e+17, 1.39e-1],
        [9.27e+17, 1.32e-1],
        [1.09e+18, 1.26e-1],
        [1.31e+18, 1.21e-1],
        [1.56e+18, 1.18e-1],
        [1.86e+18, 1.16e-1],
        [2.22e+18, 1.15e-1],
        [2.66e+18, 1.14e-1],
        [3.17e+18, 1.14e-1],
        [3.79e+18, 1.13e-1],
        [4.52e+18, 1.13e-1],
        [5.40e+18, 1.13e-1],
        [6.45e+18, 1.13e-1],
        [7.39e+18, 1.12e-1]])

    plt.plot(fiducial1a[:, 0], fiducial1a[:, 1],
             label="Model 1a (Fiducial)", c=p1[0].get_color())
    plt.plot(fiducial1b[:, 0], fiducial1b[:, 1],
             label="Model 1b (Fiducial)", c=p2[0].get_color())
    plt.legend(bbox_to_anchor=(0.4, 0.4))
    plt.show()

if __name__ == '__main__':
    main()
