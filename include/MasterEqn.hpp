#ifndef MASTEREQN_HPP
#define MASTEREQN_HPP

#include <Amplitude.hpp>

class MasterEqn {
 public:
  MasterEqn(int dim, const Amplitude *A);
  ~MasterEqn();
  void addCoupling(int m, int n, Amplitude a);
  void addDecay(int into, int outOf, double gamma);
  double getTime() const;
  void takeStep();
  const Amplitude* getState() const;

 private:
  struct MasterEqnImpl;
  MasterEqnImpl* impl;
};

#endif

