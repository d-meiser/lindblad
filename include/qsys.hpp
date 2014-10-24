#ifndef QSYS_HPP
#define QSYS_HPP

#include <complex>

typedef std::complex<double> Amplitude;

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B);
void rightApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
                Amplitude *B);

class Coupling {
 public:
  Coupling(int m, int n, Amplitude g) : m(m), n(n), g(g) {}
  void apply(int dim, const Amplitude *A, Amplitude *B);

 private:
  int m;
  int n;
  Amplitude g;
};

#endif

