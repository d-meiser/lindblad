static const char help[] =
    "Computation of steady state of a master eqaution.\n"
    "Usage:\n"
    "SteadyStateEIT:\n"
    "Compute steady state with default parameters\n"
    "SteadyStateEIT OmegaR OmegaB Delta gamma Gamma deltaB:\n"
    "Compute steady state with given parameters. All parameters must be "
    "provided\n"
    "in this form. Provide a floating point number for each parameters.";
#define _USE_MATH_DEFINES
#include <MasterEqn.hpp>
#include <Amplitude.hpp>

#include <petscmat.h>
#include <petscksp.h>

static const PetscInt N = 4;

struct SystemParameters {
  PetscReal OmegaR; /**< Rabi Freq. of Probe Beam (Hz) */
  PetscReal OmegaB; /**< Larmour Freq. of Axial Field (Hz) */
  PetscReal Delta;  /**< Detuning of Probe beam (Hz) */
  PetscReal gamma;  /**< 1/(Atom Interation Time with Probe) (Hz) */
  PetscReal Gamma;  /**< Natural Linewidth of Excited State (Hz) */
  PetscReal deltaB; /**< Larmour Freq. of Transverse Field (Hz) */
};

static PetscErrorCode getParameters(int argn, const char** argv,
                                    SystemParameters* parameters);
static PetscErrorCode readParameter(int n, const char** argv, PetscReal* param,
                                    const char* name);
static PetscErrorCode printDensityMatrix(Vec x);

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
int main(int argn, const char** argv) {
  KSP              ksp;
  PC               pc;
  Vec              x, y;
  Mat              A;
  PetscInt         i, m, n;
  PetscErrorCode   ierr;
  MasterEqn        meqn(N);
  PetscScalar      *xarr, trace;
  SystemParameters params;

  PetscFunctionBegin;
  PetscInitialize(&argn, (char***)&argv, 0, help);

  getParameters(argn, argv, &params);
  
  // intialize solution vector with initial guess for solve
  ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
  ierr = VecSetSizes(x, PETSC_DECIDE, N * N);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  for (i = 0; i < N; ++i) {
    xarr[i + i * N] = 1.0 / N;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x, &n);CHKERRQ(ierr);

  // populate y vector with source terms
  ierr = VecCreate(PETSC_COMM_WORLD, &y);CHKERRQ(ierr);
  ierr = VecSetSizes(y, PETSC_DECIDE, N * N);CHKERRQ(ierr);
  ierr = VecSetFromOptions(y);CHKERRQ(ierr);
  ierr = VecZeroEntries(y);CHKERRQ(ierr);
  ierr = VecGetLocalSize(y, &m);CHKERRQ(ierr);

  // State Labels: 0->{1,1}, 1->{1,0}, 2->{1,-1}, 3->{0,0}
  // define couplings from the Hamiltonian
  meqn.addCoupling(0, 0, params.OmegaB / 2.0);
  meqn.addCoupling(0, 1, params.deltaB / 2.0 / sqrt(2.0));
  meqn.addCoupling(1, 2, params.deltaB / 2.0 / sqrt(2.0));
  meqn.addCoupling(2, 2, -params.OmegaB / 2.0);
  meqn.addCoupling(0, 3, 1.0 / 4.0 * params.OmegaR / sqrt(6.0));
  meqn.addCoupling(2, 3, -1.0 / 4.0 * params.OmegaR / sqrt(6.0));
  meqn.addCoupling(3, 3, -params.Delta);
  // define the decays from relaxations
  // The decoherence of the excited state is probably negligible
  // compared to Gamma, but we keep it here for consistency
  meqn.addDecay(3, 3, params.gamma);
  meqn.addDecay(0, 3, params.Gamma / 3.0);
  meqn.addDecay(1, 3, params.Gamma / 3.0);
  meqn.addDecay(2, 3, params.Gamma / 3.0);
  // Mixing of ground state levels due to atoms entering and leaving
  // beam.
  meqn.addDecay(0, 0, params.gamma);
  meqn.addDecay(0, 1, params.gamma);
  meqn.addDecay(0, 2, params.gamma);
  meqn.addDecay(1, 0, params.gamma);
  meqn.addDecay(1, 1, params.gamma);
  meqn.addDecay(1, 2, params.gamma);
  meqn.addDecay(2, 0, params.gamma);
  meqn.addDecay(2, 1, params.gamma);
  meqn.addDecay(2, 2, params.gamma);

  ierr = MatCreateShell(PETSC_COMM_WORLD, m, n, N * N, N * N, &meqn, &A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))mult);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = PCSetType(pc, PCNONE);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, 1.0e-15, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
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
  for (i = 0; i < N * N; ++i) {
    xarr[i] = xarr[i] / trace;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = printDensityMatrix(x);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFinalize();
}

#undef __FUNCT__
#define __FUNCT__ "getParameters"
PetscErrorCode getParameters(int argn, const char** argv,
                             SystemParameters* parameters) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (argn == 1) {
    parameters->OmegaR = 1.25e6 * 2.0 * M_PI;
    parameters->OmegaB = 700.0e3 * 2.0 * M_PI * 0.01;
    parameters->Delta = 0.0;
    parameters->gamma = 1.0e3 * 2.0 * M_PI;
    parameters->Gamma = 6.0e6 * 2.0 * M_PI;
    /* 700kHz/Gauss, 10 milli Gauss */
    // deltaB = 700.0e3 * 2.0 * M_PI * 0.01;
    parameters->deltaB = 0.0;
  } else if (argn < 7) {
    SETERRQ(PETSC_COMM_WORLD, 1, "Insufficient parameters provided.");
  } else {
    ierr = readParameter(1, argv, &parameters->OmegaR, "OmegaR");CHKERRQ(ierr);
    ierr = readParameter(2, argv, &parameters->OmegaB, "OmegaB");CHKERRQ(ierr);
    ierr = readParameter(3, argv, &parameters->Delta, "Delta");CHKERRQ(ierr);
    ierr = readParameter(4, argv, &parameters->gamma, "gamma");CHKERRQ(ierr);
    ierr = readParameter(5, argv, &parameters->Gamma, "Gamma");CHKERRQ(ierr);
    ierr = readParameter(6, argv, &parameters->deltaB, "deltaB");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "readParameter"
PetscErrorCode readParameter(int n, const char** argv, PetscReal* param,
                    const char* name) {
  PetscFunctionBegin;
  if (sscanf(argv[n], "%lf", param) != 1) {
    SETERRQ2(PETSC_COMM_WORLD, 1, "Failed to read argument %d: %s\n", n, name);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "printDensityMatrix"
PetscErrorCode printDensityMatrix(Vec x) {
  PetscErrorCode    ierr;
  const PetscScalar *rho;
  int               i, j;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(x, &rho);CHKERRQ(ierr);
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      printf("%+1.5lf %+1.5lf  ", PetscRealPart(rho[i * N + j]), PetscImaginaryPart(rho[i * N + j]));
    }
    printf("\n");
  }
  ierr = VecRestoreArrayRead(x, &rho);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
