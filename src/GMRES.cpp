#include <GMRES.hpp>
#include <algorithm>
#include <limits>
#include <iostream>

static double absSquared(const Amplitude a) {
  return a.real() * a.real() + a.imag() * a.imag();
}

GMRES::GMRES(int dim)
    : y(m), v(dim * m), r(dim), x(dim), w(dim), bhat(m + 1),
      h((m + 1) * m), rMat((m + 1) * m), c(m), s(m) {}

void GMRES::axpy(double alpha,
                 void (*Ax)(int dim, const Amplitude *, Amplitude *, void *ctx),
                 const Amplitude *x, const Amplitude *y, Amplitude *result,
                 void *ctx) {
  Ax(r.size(), x, result, ctx);
  for (int i = 0; i < r.size(); ++i) {
    result[i] = alpha * result[i] + y[i];
  }
}

double GMRES::norm(int dim, const Amplitude *x) const {
  double nrm = 0.0;
  for (int i = 0; i < dim; ++i) {
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

bool GMRES::smallEnough(double error) const {
  if (error < 1.0e-10) {
    return true;
  } else {
    return false;
  }
}

template <typename T>
static void assignScaled(T alpha, int n, const Amplitude *x,
                             Amplitude *y) {
  for (int i = 0; i < n; ++i) {
    y[i] = alpha * x[i];
  }
}

template <typename T>
static void addAssignScaled(T alpha, int n, const Amplitude *x,
                             Amplitude *y) {
  for (int i = 0; i < n; ++i) {
    y[i] += alpha * x[i];
  }
}

void GMRES::solve(void (*A)(int, const Amplitude *, Amplitude *, void *),
                  const Amplitude *rhs, Amplitude *x0, void *ctx) {
  double rho;
  int nr;
  int dim = r.size();
  axpy(-1.0, A, x0, rhs, &r[0], ctx);
  std::copy(x0, x0 + dim, &x[0]);
  for (int j = 0; j < MAX_RESTARTS; ++j) {
    double beta = norm(dim, &r[0]);
    assignScaled(1.0 / beta, dim, &r[0], &v[0]);
    bhat[0] = beta;
    std::fill(bhat.begin() + 1, bhat.end(), 0);
    for (int i = 0; i < m; ++i) {
      A(dim, &v[i * dim], &w[0], ctx);
      for (int k = 0; k <= i; ++k) {
        h[k * m + i] = dot(dim, &v[k * dim], &w[0]);
        Amplitude hki = h[k * m + i];
        addAssignScaled(-hki, dim, &v[k * dim], &w[0]);
      }
      h[(i + 1) * m + i] = norm(dim, &w[0]);
      assignScaled(1.0 / h[(i + 1) * m + i], dim, &w[0], &v[(i + 1) * dim]);
      rMat[i] = h[i];
      for (int k = 1; k <= i; ++k) {
        Amplitude gamma = c[k - 1] * rMat[(k - 1) * m + i] +
                          std::conj(s[k - 1]) * h[k * m + i];
        rMat[k * m + i] =
            -s[k - 1] * rMat[(k - 1) * m + i] + c[k - 1] * h[k * m + i];
        rMat[(k - 1) * m + i] = gamma;
      }
      Amplitude rii = rMat[i * m + i];
      Amplitude hipi = h[(i + 1) * m + i];
      double delta = sqrt(absSquared(rii) + absSquared(hipi));
      Amplitude mu;
      Amplitude tau;
      if (absSquared(rii) < absSquared(hipi)) {
        mu = rii / hipi;
        if (std::abs(mu) < std::numeric_limits<double>::epsilon()) {
          tau = 1;
        } else {
          tau = std::conj(mu) / std::abs(mu);
        }
      } else {
        mu = hipi / rii;
        tau = mu / std::abs(mu);
      }
      c[i] = std::abs(rii) / delta;
      s[i] = std::abs(hipi) * tau / delta;
      rMat[i * m + i] = c[i] * rii + std::conj(s[i]) * hipi;
      bhat[i + 1] = -s[i] * bhat[i];
      bhat[i] = c[i] * bhat[i];
      rho = std::abs(bhat[i + 1]);
      std::cout << "iter: " << j << ", " << i << " rho == " << rho << std::endl;
      if (smallEnough(rho)) {
        nr = i;
        goto SOL;
      }
    }
    nr = m;
    y[nr] = bhat[nr] / rMat[nr * m + nr];
  SOL:
    for (int k = nr; k >= 0; --k) {
      y[k] = bhat[k];
      for (int i = k + 1; i <= nr; ++i) {
        y[k] -= rMat[k * m + i] * y[i];
      }
      y[k] /= rMat[k * m + k];
    }
    for (int i = 0; i <= nr; ++i) {
      addAssignScaled(y[i], dim, &v[i * dim], &x[0]);
    }
    if (smallEnough(rho)) {
      std::copy(x.begin(), x.end(), x0);
      return;
    }
    axpy(-1.0, A, &x[0], rhs, &r[0], ctx);
  }
}

