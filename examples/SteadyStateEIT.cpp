static const char help[] = "Computation of steady state of a master eqaution.";
#include <MasterEqn.hpp>
#include <Amplitude.hpp>

#include <petscmat.h>
#include <petscksp.h>

static const PetscInt N = 4;

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
  PetscReal      OmegaR, OmegaB, Delta, gamma, deltaB, Gamma;
  MasterEqn      meqn(N);
  PetscScalar    *xarr, trace;

  PetscInitialize(&argn, &argv, 0, help);
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

  OmegaR = 1.25e6 * 2.0 * M_PI;
  OmegaB = 700.0e3 * 2.0 * M_PI * 0.01;
  Delta = 0.0;
  gamma = 1.0e3 * 2.0 * M_PI;
  Gamma = 6.0e6 * 2.0 * M_PI;
  /* 700kHz/Gauss, 10 milli Gauss */
  //deltaB = 700.0e3 * 2.0 * M_PI * 0.01;
  deltaB = 0.0;
  meqn.addCoupling(0, 0, OmegaB / 8.0);
  meqn.addCoupling(0, 1, deltaB / 8.0 / sqrt(2.0));
  meqn.addCoupling(1, 2, deltaB / 8.0 / sqrt(2.0));
  meqn.addCoupling(2, 2, -OmegaB / 8.0);
  meqn.addCoupling(0, 3, Amplitude(1.0, -1.0) / 8.0 * OmegaR / sqrt(3.0));
  meqn.addCoupling(2, 3, -Amplitude(1.0, 1.0) / 8.0 * OmegaR / sqrt(3.0));
  meqn.addCoupling(3, 3, -Delta);
  meqn.addDecay(0, 0, gamma);
  meqn.addDecay(1, 1, gamma);
  meqn.addDecay(2, 2, gamma);
  meqn.addDecay(3, 3, gamma);
  meqn.addDecay(0, 3, Gamma);
  meqn.addDecay(1, 3, Gamma);
  meqn.addDecay(2, 3, Gamma);
  ierr = MatCreateShell(PETSC_COMM_WORLD, m, n, N * N, N * N, &meqn, &A);CHKERRQ(ierr);
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
  for (i = 0; i < N; ++i) {
    trace += xarr[i + i * N];
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "trace == %lf + i %lf\n", trace.real(), trace.imag());CHKERRQ(ierr);
  for (i = 0; i < N * N; ++i) {
    xarr[i] = xarr[i] / trace;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFinalize();
}
