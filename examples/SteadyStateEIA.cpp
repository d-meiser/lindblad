static const char help[] = "Computation of steady state of a master eqaution.";
#include <MasterEqn.hpp>
#include <Amplitude.hpp>

#include <petscmat.h>
#include <petscksp.h>

static const PetscInt N = 5 + 7 + 7;

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
  PetscInt       i, m, n, mp;
  PetscErrorCode ierr;
  PetscReal      Delta, OmegaR, OmegaB, gamma;
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
  ierr = MatCreateShell(PETSC_COMM_WORLD, m, n, N * N, N * N, &meqn, &A);CHKERRQ(ierr);

  OmegaR = 1.0;
  OmegaB = 0.00;
  Delta = -0.1;
  gamma = 1.0;
  // magnetic field interactions in the ground state manifold
  meqn.addCoupling(0, 1, OmegaB / 8.0);
  meqn.addCoupling(1, 2, OmegaB / 8.0);
  meqn.addCoupling(2, 3, -OmegaB / 8.0);
  meqn.addCoupling(3, 4, -OmegaB / 8.0);
  // magnetic field interactions in the excited state manifold
  meqn.addCoupling(5 + 0, 5 + 1, OmegaB / 8.0);
  meqn.addCoupling(5 + 1, 5 + 2, OmegaB / 8.0);
  meqn.addCoupling(5 + 2, 5 + 3, OmegaB / 8.0);
  meqn.addCoupling(5 + 3, 5 + 4, -OmegaB / 8.0);
  meqn.addCoupling(5 + 4, 5 + 5, -OmegaB / 8.0);
  meqn.addCoupling(5 + 5, 5 + 6, -OmegaB / 8.0);
  // Laser coupling (must figure out matrix elements)
  for (m = 0; m < 5; ++m) {
    for (mp = 0; mp < 3; ++mp) {
      meqn.addCoupling(m, 5 + m + mp, OmegaR);
    }
  }
  // Decay rates (must figure out branching ratios)
  for (mp = 0; mp < 7; ++mp) {
    for (m = -1; m < 2; ++m) {
      if (mp + m > 0 && mp + m < 5) {
        meqn.addDecay(m + mp, 5 + mp, gamma);
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
