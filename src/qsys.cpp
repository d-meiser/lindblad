#include <qsys.hpp>
#include <RK4.hpp>
#include <MasterEqnRhs.hpp>

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

struct MasterEqnRhsContext {
  int dim;
  const MasterEqnRhs* rhs;
};

static void applyRhs(double* x, double* y, double t, void* ctx) {
  MasterEqnRhsContext* meCtx = (MasterEqnRhsContext *)ctx;
  meCtx->rhs->apply(meCtx->dim, (const Amplitude*) x, (Amplitude*) y);
}

MasterEqn::MasterEqn(int d, const Amplitude *A)
     {
  rhs = new MasterEqnRhs();
  ctx = new MasterEqnRhsContext;
  ctx->dim = d;
  ctx->rhs = rhs;
  integrator = new RK4(2 * d * d, 0, (const double *)A, &applyRhs, 1.0e-2);
}

MasterEqn::~MasterEqn() {
  delete rhs;
  delete ctx;
  delete integrator;
}

double MasterEqn::getTime() const {
  return integrator->getTime();
}

void MasterEqn::takeStep() {
  integrator->takeStep(ctx);
}

const Amplitude* MasterEqn::getState() const {
  return (const Amplitude*) integrator->getState();
}

void MasterEqn::addCoupling(Coupling c) {
  rhs->addCoupling(c);
}

void MasterEqn::addDecay(Decay c) {
  rhs->addDecay(c);
}
