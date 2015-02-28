/*
Copyright 2014-2015 Dominic Meiser

This file is part of lindblad.

lindblad is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

lindblad is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with lindblad.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MASTEREQN_HPP
#define MASTEREQN_HPP

#include <Amplitude.hpp>
#include <SparseMatrix.hpp>
#include <LindbladExport.h>

namespace Lindblad {

/**
@brief The representation of Master Equations in Lindblad

The master equation describes the evolution of an open quantum system.
It can be written in the form

\f$
\frac{d\hat\rho(t)}{dt} = \frac{1}{i\hbar} \left[ \hat H, \hat \rho(t)\right] +
                          L\left[\hat \rho(t)\right]\;,
\f$
with \f$\hat H\f$ the Hamiltonian of the system and with \f$ L\f$
describing the dissipative parts of the dynamics.  The right hand side
of the Master equation is of the so called Lindblad form.  As a
consequence, the Master equation preserves norm and Hermiticity of the
density matrix.
*/
class LINDBLAD_API MasterEqn {
 public:
/**
@brief Construct a Master equation for a dim dimensional system.
*/
  MasterEqn(int dim);
  MasterEqn(const MasterEqn& other);
  MasterEqn& operator=(const MasterEqn& rhs);
  ~MasterEqn();

/**
@brief Add a coherent coupling between two bare states

This method adds the term
\f$
\hat H_c = (g |m\rangle\langle n| + g^* |n\rangle\langle m|) / 2
\f$
to the Hamiltonian of the system.  \f$ m = n\f$ is allowed and corresponds to a shift of an energy level by \f$ g\f$.

@param m Quantum level of first coupled state. Must be less than dim.
@param n Quantum level of second coupled state. Must be less than dim.
@param g The coupling constant.
*/
  void addCoupling(int m, int n, Amplitude g);

/**
@brief Add a decay between two bare states.

This method adds the term
\f$
-\frac{\gamma}{2}\left(
|outOf\rangle \langle outOf | \hat \rho +
\hat \rho |outOf\rangle \langle outOf | -
2|into\rangle\langle outOf|\hat \rho |outOf\rangle \langle into|\right)
\f$
To the right hand side of the Master eqaution.  The case \f$into = outOf
\f$ is allowed and corresponds to pure dephasing such as in e.g. elastic
Rayleigh scattering.

@param into  State into which the system decays.
@param outOf State out of which the system decays.
@param gamma The decay rate.

@sa addGeneralDecayOperator
*/
  void addDecay(int into, int outOf, double gamma);

/**
@brief Add a decay process between arbitrary states.

This method is a generalization of addDecay in that it allows the description of decay processes between states that are not bare energy levels.  For instance this could be used to model decay between dressed states.  Another common scenario is the case of several indistinguishable decay processes.  This method adds the term
\f$
-\frac{1}{2}\left(
\hat \lambda^\dagger\hat \lambda \hat \rho +
\hat \rho \hat \lambda^\dagger\hat \lambda -
2\hat \lambda \hat \rho \hat \lambda^\dagger
\right)
\f$
to the right hand side of the master equation.

@param lambda The decay operator.

@sa addDecay, SparseMatrix
*/
  void addGeneralDecayOperator(SparseMatrix lambda);

/**
@brief Add a particle source

This method adds a particle source for particles in a specific state
\f$into\f$.  Particles are added to the system with rate \f$gamma\f$.
To preserve the trace of the density matrix it is necessary that
particles in any state leave the system with a total rate of
\f$gamma\f$.  Sinks are always necessarily linked to sources in this
way.  Hence the slightly clumsy name of this method.  Note that a
source-sink leads to dephasing of the system at a rate \f$gamma\f$.
This dephasing process corresponds to transit time broadening.

@param into  State in which particles are added to the system
@param gamma Rate at which particles are added to the system
*/
  void addSourceSink(int into, double gamma);

/**
@brief Apply Lindblad operator to a density matrix

This evaluates the right hand side of the master equation by applying
the Lindblad super operator to a density matrix. 

@param A The density matrix to which to apply the super operator.  Must
         be an array of minimum size \f$dim^2\f$.
@param B The result.  Must be an array of minimum size \f$dim^2\f$.
*/
  void apply(const Amplitude* A, Amplitude* B) const;

/**
@brief Returns the dimension of the underlying quantum system
*/
  int getDim() const;

/**
@brief Get the level shifts of the bare levels

This returns the level shifts \f$\omega_i\f$ corresponding to the
diagonal part of the Hamiltonian of the system is \f$\sum_{i=1}^dim
\hbar\omega_i |i\rangle\langle i|\f$.

@param omegas Array with level shifts.  Must be at least of size \f$
dim\f$.
*/
  void getEnergyLevels(Amplitude* omegas) const;

 private:
  struct Impl;
  Impl* impl;
};

}

#endif

