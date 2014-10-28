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
  void apply(int dim, const Amplitude *A, Amplitude *B);

 private:
  int m;
  int n;
  Amplitude g;
};

class Decay {
 public:
  Decay(int into, int outof, double gamma)
      : into(into), outof(outof), gamma(gamma) {}
  void apply(int dim, const Amplitude *A, Amplitude *B);

 private:
  int into;
  int outof;
  double gamma;
};

class MasterEqnRhs {
 public:
  void addCoupling(Coupling c);
  void addDecay(Decay d);
  void apply(int dim, const Amplitude *A, Amplitude *B);

 private:
  std::vector<Coupling> couplings;
  std::vector<Decay> decays;
};

#endif

