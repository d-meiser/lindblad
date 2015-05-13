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
static const char help[] =
  "Computation of steady state of a master eqaution by means of a\n"
  "direct solve of a dense linear system of equations.\n\n";

#include <Lindblad.hpp>
using namespace Lindblad;
#include <vector>
#include <iostream>
#include <petscmat.h>
#include <petscksp.h>


#undef __FUNCT__
#define __FUNCT__ "main"
PetscErrorCode main(int argn, char **argv) {
  MasterEqn              meqn(2);
  double                 g = 2.0, delta = 0.1, gamma = 1.0, w = 0.2;
  Mat                    mat;
  KSP                    ksp;
  PC                     pc;
  static const PetscInt  n = 2, N = n * n;
  PetscErrorCode         ierr;
  std::vector<Amplitude> matrix(N * N);
  PetscInt               i, j, m = N, idxm[N];
  Vec                    rhs, x;
  PetscScalar            *xarr;
  PetscReal              trace;
  PetscBool              print_matrix, flg;

  PetscFunctionBegin;
  PetscInitialize(&argn, &argv, 0, help);

  ierr = PetscOptionsGetReal("", "-g", &g, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal("", "-delta", &delta, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal("", "-gamma", &gamma, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal("", "-w", &w, &flg);CHKERRQ(ierr);

  meqn.addCoupling(0, 1, g);
  meqn.addCoupling(0, 0, -0.5 * delta);
  meqn.addCoupling(1, 1, 0.5 * delta);
  meqn.addDecay(0, 1, gamma);
  meqn.addDecay(1, 0, w);

  meqn.buildMatrix(&matrix[0]);
  ierr = PetscOptionsGetBool("", "-print_matrix", &print_matrix, &flg);CHKERRQ(ierr);
  if (print_matrix) {
    for (i = 0; i < N; ++i) {
      for (j = 0; j < N; ++j) {
        std::cout << matrix[i * N + j] << " ";
      }
      std::cout << std::endl;
    }
  }

  ierr = MatCreate(PETSC_COMM_WORLD, &mat);CHKERRQ(ierr);
  ierr = MatSetType(mat, MATSEQBAIJ);CHKERRQ(ierr);
  for (i = 0; i < N; ++i) {
    idxm[i] = i;
  }
  ierr = MatSetFromOptions(mat);CHKERRQ(ierr);
  ierr = MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, N, N);CHKERRQ(ierr);
  ierr = MatSetUp(mat);CHKERRQ(ierr);
  ierr = MatSetValues(mat, m, idxm, m, idxm, &matrix[0], INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = PCSetType(pc, PCILU);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, mat, mat);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecSetSizes(x, PETSC_DECIDE, N);CHKERRQ(ierr);
  ierr = VecDuplicate(x, &rhs);CHKERRQ(ierr);
  ierr = VecZeroEntries(x);CHKERRQ(ierr);
  ierr = VecZeroEntries(rhs);CHKERRQ(ierr);

  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
    xarr[i + i * n] = 1.0 / n;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);
  ierr = KSPSolve(ksp, rhs, x);CHKERRQ(ierr);

  ierr = VecGetArray(x, &xarr);CHKERRQ(ierr);
  trace = 0;
  for (i = 0; i < n; ++i) {
      trace += PetscRealPart(xarr[i * n + i]);
  }
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      xarr[i * n + j] /= trace;
    }
  }
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      std::cout << xarr[i * n + j] << " ";
    }
    std::cout << std::endl;
  }
  ierr = VecRestoreArray(x, &xarr);CHKERRQ(ierr);

  ierr = MatDestroy(&mat);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  PetscFinalize();
  PetscFunctionReturn(0);
}
