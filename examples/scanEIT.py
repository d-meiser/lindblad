import subprocess
import numpy as np
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.pyplot as plt

def get_steady_state(arguments = None):
    if arguments:
        output = subprocess.check_output(["./SteadyStateEIT"] +
                [str(a) for a in arguments])
    else:
        output =  subprocess.check_output(["./SteadyStateEIT"])
    density_matrix = np.matrix(output)
    density_matrix = density_matrix.reshape(4, 8)
    density_matrix = density_matrix[:, 0:8:2] + 1.0j * density_matrix[:,1:8:2]
    return density_matrix

def get_polarization(arguments = None):
    density_matrix = get_steady_state(arguments)
    return density_matrix[0, 3]

OmegaR = 1.25e6 * 2.0 * np.pi
Delta = 0.0;
gamma = 1.0e3 * 2.0 * np.pi
Gamma = 6.0e6 * 2.0 * np.pi
deltaB = np.arange(-0.1, 0.1, 0.0005) * 700.0e3 * 2.0 * np.pi

OmegaB = 700.0e3 * 2.0 * np.pi * 0.01
absorptionsNonZeroField = np.imag(np.array([
  get_polarization([OmegaR, OmegaB, Delta, gamma, Gamma, db])
                    for db in deltaB]))

OmegaB = 700.0e3 * 2.0 * np.pi * 0.0
absorptionsZeroField = np.imag(np.array([
  get_polarization([OmegaR, OmegaB, Delta, gamma, Gamma, db])
                    for db in deltaB]))

plt.plot(deltaB / (700.0e3 * 2.0 * np.pi), 1.0e3 * absorptionsZeroField)
plt.plot(deltaB / (700.0e3 * 2.0 * np.pi), 1.0e3 *
        absorptionsNonZeroField,'--r')

plt.xlabel(r'$B_z({\rm G})$')
plt.ylabel(r'${\rm Absorption\; [arb. units]}$')
plt.gcf().set_size_inches(4, 3)
plt.gcf().subplots_adjust(bottom = 0.15, left = 0.15, top = 0.97, right = 0.95)
plt.savefig('absorption.pdf')

