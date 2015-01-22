#!/usr/bin/python
import sys
import subprocess
import numpy as np

def get_evolution(dim, arguments = None):
    """Returns density matrix evolution for given parameters.

    The evolution is computed using the EIT program.

    dim       Dimension of the density matrix.
    arguments List of arguments for to be passed to EIT.
    """

    if not arguments:
        raise RuntimeError("Must provide arguments")
    if not len(arguments) == 9:
        raise RuntimeError("Wrong number of arguments")

    output = subprocess.check_output(["./EIT"] +
                [str(a) for a in arguments])
    density_matrix_evolution = np.matrix(output)
    density_matrix_evolution = density_matrix_evolution.reshape(
            density_matrix_evolution.shape[1] / (2 * dim**2),
            2 * dim**2)
    density_matrix_evolution = density_matrix_evolution[:, 0:(2*dim**2):2] + 1.0j * density_matrix_evolution[:,1:(2*dim**2):2]
    return density_matrix_evolution


def get_steady_state(arguments = None):
    """Compute steady state of the density matrix for given parameters

    The number of time steps specified in the last entry in the
    arguments is taken as an initial guess.  The time interval is
    doubled started from there until the density matrix doesn't change
    much any more.

    arguments List of arguments to be passed to EIT without numDumps.
              numDumps is compute in this function.
    """

    if not arguments:
        raise RuntimeError("Must provide arguments to get_steady_state()")
    if not len(arguments) == 8:
        raise RuntimeError("Wrong number of arguments")

    # We are dealing with a 4 state system.
    dim = 4

    # The variation is a measure of how much the density matrix is still
    # changing as a function of time. We initialize the variation to a
    # "large" value to force at least one iteration.
    variation = 1.0

    # Convergence to steady state has occured when the variation is
    # smaller than tol.
    tol = 1.0e-4
    num_iter = 0
    max_num_iters = 20

    # Iterate to convergence by integrating over increasingly long time
    # intervals.
    while variation > tol and num_iter < max_num_iters:
        # Double the number of time steps
        arguments[-1] = arguments[-1] * 2

        # Obtain snap shots of the density matrix at 4 times
        num_snapshots = 4
        num_dump = arguments[-1] / num_snapshots
        density_matrix_evolution = get_evolution(dim, arguments + [num_dump])

        # Update variation by comparing density matrix at end of time
        # interval with density matrix at second to last time
        variation = np.linalg.norm(density_matrix_evolution[-1] -
                                   density_matrix_evolution[-2])

        print num_iter, arguments[-1], variation
        num_iter = num_iter + 1

    if variation > tol:
        raise RuntimeError("Failed to converge to steady state in "
                           + str(max_num_iters) + "iterations.")

    # Converged. Post process solution and return.
    steady_state = density_matrix_evolution[-1]
    steady_state = steady_state.reshape(dim, dim)
    return steady_state

def main(argv):
    OmegaR = 1.25e6 * 2.0 * np.pi
    Delta = 0.0;
    gamma = 1.0e3 * 2.0 * np.pi
    Gamma = 6.0e6 * 2.0 * np.pi
    OmegaB = np.arange(-0.1, 0.1, 0.0005) * 700.0e3 * 2.0 * np.pi

    deltaB = 700.0e3 * 2.0 * np.pi * 0.01

    dt = 4.0e-8
    num_steps_guess = 10000
    steady_state = get_steady_state([OmegaR, OmegaB[0], Delta, gamma, Gamma,
            deltaB, dt, num_steps_guess])
    print steady_state

if __name__== "__main__":
    main(sys.argv)

