#ifndef QSYS_HPP
#define QSYS_HPP

#include <complex>

typedef std::complex<double> Amplitude;

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B);

void rightApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B);
#endif

