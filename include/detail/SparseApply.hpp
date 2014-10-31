#ifndef SPARSE_APPLY_HPP
#define SPARSE_APPLY_HPP

#include <qsys.hpp>

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B);
void rightApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
                Amplitude *B);

#endif
