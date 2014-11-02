#ifndef MASTER_EQN_EVOLUTION_HPP
#define MASTER_EQN_EVOLUTION_HPP

#include <MasterEqn.hpp>

class MasterEqnEvolution {
 public:
  MasterEqnEvolution(const MasterEqn& eqn, const Amplitude* initialState);
  ~MasterEqnEvolution();
  double getTime() const;
  void takeStep();
  const Amplitude* getState() const;

  friend struct MasterEqnEvolutionContext;

 private:
  struct Impl;
  Impl* impl;
};

#endif

