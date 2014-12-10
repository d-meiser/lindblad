#include <GMRES.hpp>

GMRES::GMRES(int dim) : y(dim * m), r(dim) {}

void GMRES::axpy(double alpha,
                 void (*Ax)(int dim, const Amplitude *, Amplitude *, void *ctx),
                 const Amplitude *x, const Amplitude *y, Amplitude *result,
                 void *ctx) {
  Ax(r.size(), x, result, ctx);
  for (int i = 0; i < r.size(); ++i) {
    result[i] = alpha * result[i] + y[i];
  }
}

void GMRES::solve(void (*Ax)(int dim, const Amplitude *x, Amplitude *y,
                             void *ctx),
                  const Amplitude *rhs, Amplitude *x) {}

