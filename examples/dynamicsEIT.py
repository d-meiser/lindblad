#!/usr/bin/python
import sys
import subprocess
import numpy as np
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.pyplot as plt

def get_dynamics(arguments = None):
    if arguments:
        output = subprocess.check_output(["./EIT"] +
                [str(a) for a in arguments])
    else:
        output =  subprocess.check_output(["./EIT"])
    density_matrix = np.matrix(output)
    print density_matrix
    density_matrix = density_matrix.reshape(density_matrix.shape[1] / 8, 8)
    density_matrix = density_matrix[:, 0:8:2] + 1.0j * density_matrix[:,1:8:2]
    return density_matrix

def main(argv):
    OmegaR = 1.25e6 * 2.0 * np.pi
    Delta = 0.0;
    gamma = 1.0e3 * 2.0 * np.pi
    Gamma = 6.0e6 * 2.0 * np.pi
    OmegaB = np.arange(-0.1, 0.1, 0.0005) * 700.0e3 * 2.0 * np.pi

    deltaB = 700.0e3 * 2.0 * np.pi * 0.01

    num_steps = 10000
    num_dump = 100
    num_snapshots = num_steps / num_dump
    dt = 4.0e-8
    dynamics = get_dynamics([OmegaR, OmegaB[0], Delta, gamma, Gamma,
            deltaB, dt, num_steps, num_dump])
    populations = [np.array([dynamics[i*4+n,n].real for i in
                            range(num_snapshots)]) for n in range(4)]
    plt.clf()
    plots=[]
    for n in range(4):
        plots.append(plt.plot(populations[n])[0])
    plt.legend(plots,
            [r'$\rho_{0,0}$', r'$\rho_{1,1}$', r'$\rho_{2,2}$', r'$\rho_{3,3}$'])
    plt.xlabel(r'$t/\mu s$')
    plt.ylabel(r'$\rho_{n,n}(t)$')
    plt.gcf().set_size_inches(4, 3)
    plt.gcf().subplots_adjust(bottom = 0.15, left = 0.2, top = 0.97, right = 0.95)
    plt.savefig("populations.png")

if __name__== "__main__":
    main(sys.argv)
