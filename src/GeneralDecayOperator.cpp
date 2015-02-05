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
#include <GeneralDecayOperator.hpp>
#include <SparseApply.hpp>
#include <Memory.hpp>

GeneralDecayOperator::GeneralDecayOperator(SparseMatrix lambda)
    : lambda(lambda) {
  for (std::vector<SparseMatrixEntry>::const_iterator e =
           lambda.entries.begin();
       e != lambda.entries.end(); ++e) {
    lambdaDagger.add(SparseMatrixEntry(e->col, e->row, std::conj(e->element)));
  }
  for (std::vector<SparseMatrixEntry>::const_iterator eDagger =
           lambdaDagger.entries.begin();
       eDagger != lambdaDagger.entries.end(); ++eDagger) {
    for (std::vector<SparseMatrixEntry>::const_iterator e =
             lambda.entries.begin();
         e != lambda.entries.end(); ++e) {
      if (eDagger->col == e->row) {
        lambdaDaggerLambda.add(SparseMatrixEntry(
            eDagger->row, e->col, eDagger->element * e->element));
      }
    }
  }
}

void GeneralDecayOperator::apply(int dim, const Amplitude *A,
                                 Amplitude *B) const {
  for (std::vector<SparseMatrixEntry>::const_iterator e =
           lambdaDaggerLambda.entries.begin();
       e != lambdaDaggerLambda.entries.end(); ++e) {
    leftApply(e->row, e->col, -0.5 * e->element, dim, A, B);
    rightApply(e->row, e->col, -0.5 * e->element, dim, A, B);
  }
  std::vector<Amplitude> tmp(dim * dim, 0);
  for (std::vector<SparseMatrixEntry>::const_iterator e =
           lambda.entries.begin();
       e != lambda.entries.end(); ++e) {
    leftApply(e->row, e->col, e->element, dim, A, &tmp[0]);
  }
  for (std::vector<SparseMatrixEntry>::const_iterator e =
           lambdaDagger.entries.begin();
       e != lambdaDagger.entries.end(); ++e) {
    rightApply(e->row, e->col, e->element, dim, &tmp[0], B);
  }
}

