#ifndef COUPLING_HPP
#define COUPLING_HPP

#include <Amplitude.hpp>

class Coupling {
 public:
  Coupling(int m, int n, Amplitude g) : m(m), n(n), g(g) {}
  void apply(int dim, const Amplitude *A, Amplitude *B) const;

 private:
  int m;
  int n;
  Amplitude g;
};

#endif

