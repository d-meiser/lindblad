#include <qsys.hpp>

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B) {
  B += row * dim;
  A += row * dim;
  for (int c = 0; c < dim; ++c) {
    B[c] += alpha * A[c];
  }
}
