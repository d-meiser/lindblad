static const char help[] = "Computation of steady state of a master eqaution.";
#include <MasterEqn.hpp>
#include <Amplitude.hpp>

using namespace Lindblad;

#include <petscmat.h>
#include <petscksp.h>

static const PetscInt numAtoms = 4;
static const PetscInt numPhotons = 10;
static const PetscInt dim = (numPhotons + 1) * (1 << numAtoms);

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
  PetscInt       i, m, n, mp, mpp;
  PetscErrorCode ierr;
  PetscReal      w, gamma, g, kappa;
  MasterEqn      meqn(dim);
  PetscScalar    *xarr, trace;

  PetscInitialize(&argn, &argv, 0, help);
  ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
  ierr = VecSetSizes(x, PETSC_DECIDE, dim * dim);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  for (i = 0; i < dim; ++i) {
    xarr[i + i * dim] = 1.0 / dim;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &y);CHKERRQ(ierr);
  ierr = VecSetSizes(y, PETSC_DECIDE, dim * dim);CHKERRQ(ierr);
  ierr = VecSetFromOptions(y);CHKERRQ(ierr);
  ierr = VecZeroEntries(y);CHKERRQ(ierr);
  ierr = VecGetLocalSize(y, &m);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x, &n);CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD, m, n, dim * dim, dim * dim, &meqn, &A);CHKERRQ(ierr);

  gamma = 1.0;
  w = 5.0 * gamma;
  g = 1.0;
  kappa = 10.0;
  for (n = 0; n < numPhotons + 1; ++n) {
    for (m = 0; m < numAtoms; ++m) {
      for (mp = 0; mp < 1 << (m - 1); ++mp) {
        for (mpp = 1 << m; mpp < 1 << numAtoms; ++mpp) {
          meqn.addDecay(n * (1 << numAtoms) + (mp | mpp),
                        n * (1 << numAtoms) + (mp | (1 << m) | mpp), gamma);
          meqn.addDecay(n * (1 << numAtoms) + (mp | (1 << m) | mpp),
                        n * (1 << numAtoms) + (mp | mpp), w);
        }
      }
    }
  }
  for (n = 0; n < numPhotons; ++n) {
    for (m = 0; m < 1 << numAtoms; ++m) {
      meqn.addDecay(n * (1 << numAtoms) + m, (n + 1) * (1 << numAtoms) + m,
                    kappa);
    }
  }
  for (n = 0; n < numPhotons; ++n) {
    for (m = 0; m < numAtoms; ++m) {
      for (mp = 0; mp < 1 << (m - 1); ++mp) {
        for (mpp = 1 << m; mpp < 1 << numAtoms; ++mpp) {
          meqn.addCoupling((n + 1) * (1 << numAtoms) + (mp | mpp),
                           n * (1 << numAtoms) + (mp | (1 << m) | mpp),
                           std::sqrt(n + 1) * g);
        }
      }
    }
  }

  ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))mult);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = PCSetType(pc, PCNONE);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSolve(ksp, y, x);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  trace = 0;
  for (i = 0; i < dim; ++i) {
    trace += xarr[i + i * dim];
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "trace == %lf + i %lf\n", trace.real(), trace.imag());CHKERRQ(ierr);
  for (i = 0; i < dim * dim; ++i) {
    xarr[i] = xarr[i] / trace;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFinalize();
}
