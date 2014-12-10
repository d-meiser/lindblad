#include <GMRES.hpp>
#include <algorithm>

GMRES::GMRES(int dim) : y(dim * m), r(dim), x(dim) {}

void GMRES::axpy(double alpha,
                 void (*Ax)(int dim, const Amplitude *, Amplitude *, void *ctx),
                 const Amplitude *x, const Amplitude *y, Amplitude *result,
                 void *ctx) {
  Ax(r.size(), x, result, ctx);
  for (int i = 0; i < r.size(); ++i) {
    result[i] = alpha * result[i] + y[i];
  }
}

double GMRES::norm(const std::vector<Amplitude>& x) {
  double nrm = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    nrm += x[i].real() * x[i].real() + x[i].imag() * x[i].imag();
  }
  return sqrt(nrm);
}

void GMRES::solve(void (*A)(int dim, const Amplitude *, Amplitude *,
                             void *ctx),
                  const Amplitude *rhs, Amplitude *x0,  void *ctx) {
  axpy(-1.0, A, x0, rhs, &r[0], ctx);
  std::copy(x0, x0 + r.size(), &x[0]);
  for (int j = 0; j < MAX_RESTARTS; ++j) {
    double beta = norm(r);
    for (int jj = 0; j < r.size(); ++j) {
      y[jj] = r[jj] / beta;
    }
  }
}

