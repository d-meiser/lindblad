#ifndef QSYS_HPP
#define QSYS_HPP

#include <complex>
#include <vector>

typedef std::complex<double> Amplitude;

void leftApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
               Amplitude *B);
void rightApply(int row, int col, Amplitude alpha, int dim, const Amplitude *A,
                Amplitude *B);

class Coupling {
 public:
  Coupling(int m, int n, Amplitude g) : m(m), n(n), g(g) {}
  void apply(int dim, const Amplitude *A, Amplitude *B) const;

 private:
  int m;
  int n;
  Amplitude g;
};

class Decay {
 public:
  Decay(int into, int outof, double gamma)
      : into(into), outof(outof), gamma(gamma) {}
  void apply(int dim, const Amplitude *A, Amplitude *B) const;

 private:
  int into;
  int outof;
  double gamma;
};

class MasterEqnRhs {
 public:
  void addCoupling(Coupling c);
  void addDecay(Decay d);
  void apply(int dim, const Amplitude *A, Amplitude *B) const;

 private:
  std::vector<Coupling> couplings;
  std::vector<Decay> decays;
};

class Integrator;
struct MasterEqnRhsContext;
class MasterEquation {
 public:
  MasterEquation(int dim, const Amplitude *A, const MasterEqnRhs *rhs);
  ~MasterEquation();
  double getTime() const;
  void takeStep();
  const Amplitude* getState() const;

 private:
  MasterEqnRhsContext *ctx;
  Integrator* integrator;
};

#endif

