#ifndef DECAY_HPP
#define DECAY_HPP

#include <Amplitude.hpp>

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

#endif

