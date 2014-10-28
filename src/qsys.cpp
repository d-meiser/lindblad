#include <qsys.hpp>

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B) {
  B += row * dim;
  A += col * dim;
  for (int c = 0; c < dim; ++c) {
    B[c] += alpha * A[c];
  }
}

void rightApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
                Amplitude *B) {
  B += col;
  A += row;
  for (int r = 0; r < dim; ++r) {
    B[r * dim] += alpha * A[r * dim];
  }
}

void Coupling::apply(int dim, const Amplitude *A, Amplitude *B) const {
  static const Amplitude I(0, 1);
  leftApply(m, n, 0.5 * g / I, dim, A, B);
  leftApply(n, m, 0.5 * conj(g) / I, dim, A, B);
  rightApply(m, n, -0.5 * g / I, dim, A, B);
  rightApply(n, m, -0.5 * conj(g) / I, dim, A, B);
}

void Decay::apply(int dim, const Amplitude *A, Amplitude *B) const {
  leftApply(outof, outof, -0.5 * gamma, dim, A, B);
  rightApply(outof, outof, -0.5 * gamma, dim, A, B);
  B[into + into * dim] += gamma * A[outof + outof * dim];
}

void MasterEqnRhs::addCoupling(Coupling c) {
  couplings.push_back(c);
}

void MasterEqnRhs::addDecay(Decay d) {
  decays.push_back(d);
}

void MasterEqnRhs::apply(int dim, const Amplitude *A, Amplitude *B) const {
  for (std::vector<Coupling>::const_iterator c = couplings.begin();
       c != couplings.end(); ++c) {
    c->apply(dim, A, B);
  }
  for (std::vector<Decay>::const_iterator d = decays.begin(); d != decays.end();
       ++d) {
    d->apply(dim, A, B);
  }
}

MasterEquation::MasterEquation(int d, const Amplitude *A, const MasterEqnRhs *r)
    : dim(d), time(0), dt(1.0e-3), state(A, A + dim * dim), rhs(r) {}

void MasterEquation::takeStep() { time += dt; }
