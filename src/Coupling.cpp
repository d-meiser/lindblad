#include <Coupling.hpp>
#include <SparseApply.hpp>

void Coupling::apply(int dim, const Amplitude *A, Amplitude *B) const {
  static const Amplitude I(0, 1);
  leftApply(m, n, 0.5 * g / I, dim, A, B);
  leftApply(n, m, 0.5 * conj(g) / I, dim, A, B);
  rightApply(m, n, -0.5 * g / I, dim, A, B);
  rightApply(n, m, -0.5 * conj(g) / I, dim, A, B);
}


