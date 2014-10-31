#ifndef MASTER_EQN_RHS_HPP
#define MASTER_EQN_RHS_HPP

#include <qsys.hpp>

class MasterEqnRhs {
 public:
  void addCoupling(Coupling c);
  void addDecay(Decay d);
  void apply(int dim, const Amplitude *A, Amplitude *B) const;

 private:
  std::vector<Coupling> couplings;
  std::vector<Decay> decays;
};

#endif

