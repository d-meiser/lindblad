/*
Copyright 2014-2015 Dominic Meiser

This file is part of lindblad.

lindblad is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

lindblad is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with lindblad.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
\par
This example illustrates how to use the Lindblad library to compute
steady state solutions of master equations.  To keep things simple we
consider the case of a driven two level atom.

\par
In this example we use the PETSc library to solve the linear system \f$
d\hat\rho / dt = 0 \f$.  PETSc has to be installed separately, see for
example the installation instructions in the file README.md.  Note that
the iterative (matrix free) methods do not always converge for the
solution of the steady state equation because good preconditioners are
difficult to construct in the general case.  This can lead to solver
failure for ill-conditioned systems with e.g. largely different time
scales.  Furthermore, this method cannot be used when the Hamiltonian is
explicitly time dependent.  While it is possible to transform time
dependencies away in many cases of practical interest this cannot be
done in general.

\par
Many aspects of this example can be configured from the command line
using PETSc options.  For instance, to get a log of the convergence of
the GMRES solver one may use

    ./SteadyState -ksp_monitor

\par
To see more options use <tt>-help</tt>.

\par
Alternative approaches are:

- Propagate an initial state to steady state using MasterEqnEvolution.
This is the most generally applicable approach.  Typically it is also
the most robust.  This method is less efficient than an iterative
solution using e.g. GMRES.  In particular, a large number of time steps
may be necessary when there are large disparaties of relaxation rates.

- Form the matrix corresponding to the Master Equation and use a direct
solver (sparse or dense).  This method may require large amounts of
storage due to the size of the operators in Liouville space.

\par
We need the headers <tt>MasterEqn.hpp</tt> and <tt>Amplitude.hpp</tt>.
*/
static const char help[] = "Computation of steady state of a simple master equation.";
#include <MasterEqn.hpp>
using namespace Lindblad;
#include <Amplitude.hpp>

/*
\par
From the PETSc library we need <tt>petscmat.h</tt> for <tt>MatShell</tt>
and <tt>petscksp.h</tt> for linear solvers:
*/
#include <petscmat.h>
#include <petscksp.h>

/*
\par
We are dealing with a two level system, i.e. the Hilbert space dimension is 2:
*/
static const PetscInt N = 2;

/*
\par
The function <tt>mult</tt> implements the application of the linear
operator corresponding to the right hand side of the master equation.
To interface more easily with PETSc, the operator is wrapped in a
<tt>MatShell</tt>.  The context of the <tt>MatShell</tt> is simply a
pointer to a <tt>MasterEqn</tt> object describing our system.

\remark Note that we are following the style conventions of PETSc in
this example.  This includes defining the <tt>__FUNCT__</tt> macro,
usage of <tt>PetscFunctionBegin</tt>, and <tt>PetscFunctionReturn</tt>,
error handling, and code formatting.
*/
#undef __FUNCT__
#define __FUNCT__ "mult"
int mult(Mat A, Vec x, Vec y) {
  PetscErrorCode     ierr;
  MasterEqn*         meqn;
  PetscScalar*       yarr;
  const PetscScalar* xarr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A, (void*)&meqn);CHKERRQ(ierr);
  ierr = VecGetArrayRead(x, &xarr);CHKERRQ(ierr);
  ierr = VecGetArray(y, &yarr);CHKERRQ(ierr);
  meqn->apply((const Amplitude*)xarr, (Amplitude*)yarr);
  ierr = VecRestoreArrayRead(x, &xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(y, &yarr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argn, char **argv) {
  KSP            ksp;
  PC             pc;
  Vec            x, y;
  Mat            A;
  PetscInt       i, m, n;
  PetscErrorCode ierr;
  PetscReal      g, delta, gamma, trace;
  MasterEqn      meqn(2);
  PetscScalar    *xarr;

  PetscInitialize(&argn, &argv, 0, help);
/*
First we create vectors for the solution, <tt>x</tt>, and right hand
side, <tt>y</tt>.  The right hand side is set to zero because we are
looking for a steady state.  The solution is initiallized to an
incoherent mixture with equal probability to find the partical in each
bare state.
*/
  ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
  ierr = VecSetSizes(x, PETSC_DECIDE, N * N);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  for (i = 0; i < N; ++i) {
    xarr[i + i * N] = 1.0 / N;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &y);CHKERRQ(ierr);
  ierr = VecSetSizes(y, PETSC_DECIDE, N * N);CHKERRQ(ierr);
  ierr = VecSetFromOptions(y);CHKERRQ(ierr);
  ierr = VecZeroEntries(y);CHKERRQ(ierr);
  ierr = VecGetLocalSize(y, &m);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x, &n);CHKERRQ(ierr);

/*
Next we set up the <tt>MasterEqn</tt> object with coupling, detuning,
and decay out of the excited state.  We also create a <tt>MatShell</tt>
and set the <tt>mult</tt> function for application of the linear
operator.
*/
  g = 2.0;
  meqn.addCoupling(0, 1, g);
  delta = 0.1;
  meqn.addCoupling(1, 1, delta);
  gamma = 1.0;
  meqn.addDecay(0, 1, gamma);
  ierr = MatCreateShell(PETSC_COMM_WORLD, m, n, N * N, N * N, &meqn, &A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))mult);CHKERRQ(ierr);

/*
Then we set up the iterative solver.  We set the preconditioner type to
<tt>PCNONE</tt> because we don't have an assembled matrix and hence
cannot use preconditioners like ILU.  We tell the KSP that the density
matrix <tt>x</tt> created above should be used as initial guess.  Note
that more options can be set using the command line.
*/
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = PCSetType(pc, PCNONE);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSolve(ksp, y, x);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

/*
We need to ensure that the trace of the steady state density matrix is 1.
*/
  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  trace = 0;
  for (i = 0; i < N; ++i) {
    trace += xarr[i + i * N].real();
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "trace == %lf\n", trace);CHKERRQ(ierr);
  for (i = 0; i < N * N; ++i) {
    xarr[i] /= trace;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

/*
Finally we clean up all remaining resources.
*/
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFinalize();
}
