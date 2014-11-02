#ifndef MASTEREQN_HPP
#define MASTEREQN_HPP

#include <Amplitude.hpp>

class MasterEqn {
 public:
  MasterEqn(int dim, const Amplitude *A);
  MasterEqn(const MasterEqn& other);
  MasterEqn& operator=(const MasterEqn& rhs);
  ~MasterEqn();
  void addCoupling(int m, int n, Amplitude a);
  void addDecay(int into, int outOf, double gamma);
  double getTime() const;
  void takeStep();
  const Amplitude* getState() const;
  void apply(int dim, const Amplitude* A, Amplitude* B) const;
  int getDim() const;

  friend struct MasterEqnRhsContext;
 private:
  struct Impl;
  Impl* impl;
};

#endif

