#include <GMRES.hpp>
#include <algorithm>

GMRES::GMRES(int dim)
    : y(dim * m), v(dim * m), r(dim), x(dim), w(dim), bhat(dim), h(m * m) {}

void GMRES::axpy(double alpha,
                 void (*Ax)(int dim, const Amplitude *, Amplitude *, void *ctx),
                 const Amplitude *x, const Amplitude *y, Amplitude *result,
                 void *ctx) {
  Ax(r.size(), x, result, ctx);
  for (int i = 0; i < r.size(); ++i) {
    result[i] = alpha * result[i] + y[i];
  }
}

double GMRES::norm(const std::vector<Amplitude>& x) const {
  double nrm = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    nrm += x[i].real() * x[i].real() + x[i].imag() * x[i].imag();
  }
  return sqrt(nrm);
}

Amplitude GMRES::dot(int dim, const Amplitude *x, const Amplitude *y) const {
  Amplitude result = 0;
  for (int i = 0; i < dim; ++i) {
    result += std::conj(x[i]) * y[i];
  }
  return result;
}

void GMRES::solve(void (*A)(int, const Amplitude *, Amplitude *, void *),
                  const Amplitude *rhs, Amplitude *x0, void *ctx) {
  int dim = r.size();
  axpy(-1.0, A, x0, rhs, &r[0], ctx);
  std::copy(x0, x0 + dim, &x[0]);
  for (int j = 0; j < MAX_RESTARTS; ++j) {
    double beta = norm(r);
    for (int jj = 0; jj < dim; ++jj) {
      v[jj] = r[jj] / beta;
    }
    bhat[0] = beta;
    std::fill(bhat.begin() + 1, bhat.end(), 0);
    for (int i = 0; i < m; ++i) {
      A(dim, &v[i * dim], &w[0], ctx);
      for (int k = 0; k < i; ++k) {
        h[k * m + i] = dot(dim, &v[k * dim], &w[0]);
        for (int jj = 0; jj < dim; ++jj) {
          w[jj] -= h[k * m + i] * v[k * dim + jj];
        }
      }
    }
  }
}

