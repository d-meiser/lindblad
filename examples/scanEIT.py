#!/usr/bin/python

import sys
import subprocess
import numpy as np
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
try:
    import joblib
    n_jobs = 4
    print "Using joblib with", n_jobs, "processes."
    have_joblib = True
except ImportError:
    print "joblib not available."
    have_joblib = False
    

if (len(sys.argv) == 2):
	# one argument has been provided for output plot name
	outputFilename = str(sys.argv[1])
elif (len(sys.argv) == 1):
	# no filename provided, use default
	outputFilename = 'test.png'
else:
	# improper usage, use default and inform
	print("Usage: python scanEIT-SmallRange.py filename")
	outputFilename = 'test.png'

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
    return [density_matrix[0, 3] , density_matrix[2, 3]]

def absorption(polarizationArray):
	# post-process density matrix elements 
	# to get absorption from atomic polarization
	return polarizationArray[:,0].imag - polarizationArray[:,1].imag
	
def rotation(polarizationArray):
	# post-process density matrix elements 
	# to get light polarization rotation from atomic polarization
	return polarizationArray[:,0].real + polarizationArray[:,1].real

def compute_polarization(OmegaR, OmegaBs, Delta, gamma, Gamma, deltaB):
    if have_joblib:
        polarizations = joblib.Parallel(n_jobs = n_jobs)(
                joblib.delayed(get_polarization)(
                    [OmegaR, ob, Delta, gamma, Gamma, deltaB])
                for ob in OmegaBs)
    else:
        polarizations = [get_polarization([OmegaR, ob, Delta, gamma, Gamma, deltaB]) for ob in OmegaBs]
    polarizations = np.array(polarizations)
    return polarizations

def main(argv):
    OmegaR = 1.25e6 * 2.0 * np.pi
    Delta = 0.0;
    gamma = 1.0e3 * 2.0 * np.pi
    Gamma = 6.0e6 * 2.0 * np.pi
    OmegaB = np.arange(-0.1, 0.1, 0.0005) * 700.0e3 * 2.0 * np.pi

    deltaB = 700.0e3 * 2.0 * np.pi * 0.01
    NonZeroPolarization = compute_polarization(OmegaR, OmegaB, Delta, gamma, Gamma, deltaB)
    absorptionNonZeroField = absorption(NonZeroPolarization)
    rotationNonZeroField = rotation(NonZeroPolarization)

    deltaB = 700.0e3 * 2.0 * np.pi * 0.0
    ZeroPolarization = compute_polarization(OmegaR, OmegaB, Delta, gamma, Gamma, deltaB)
    absorptionZeroField = absorption(ZeroPolarization)
    rotationZeroField = rotation(ZeroPolarization)

    plt.subplot(2,1,1)
    plt.plot(OmegaB / (700.0e3 * 2.0 * np.pi), 1.0e3 * absorptionZeroField)
    plt.plot(OmegaB / (700.0e3 * 2.0 * np.pi), 1.0e3 *
            absorptionNonZeroField,'--r')

    plt.ylabel(r'${\rm Absorption\; [arb. units]}$')

    plt.subplot(2,1,2)
    plt.plot(OmegaB / (700.0e3 * 2.0 * np.pi), 1.0e3 * rotationZeroField)
    plt.plot(OmegaB / (700.0e3 * 2.0 * np.pi), 1.0e3 *
            rotationNonZeroField,'--r')
            
    plt.xlabel(r'$B_z({\rm G})$')
    plt.ylabel(r'${\rm Faraday Rotation\; [arb. units]}$')

    plt.gcf().set_size_inches(4, 6)
    plt.gcf().subplots_adjust(bottom = 0.1, left = 0.2, top = 0.97, right = 0.95)
    plt.savefig(outputFilename,format='png')

if __name__== "__main__":
    main(sys.argv)
