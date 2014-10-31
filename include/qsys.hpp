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

class MasterEqn {
 public:
  MasterEqn(int dim, const Amplitude *A);
  ~MasterEqn();
  void addCoupling(Coupling c);
  void addDecay(Decay d);
  double getTime() const;
  void takeStep();
  const Amplitude* getState() const;

 private:
  struct MasterEqnImpl;
  MasterEqnImpl* impl;
};

#endif

