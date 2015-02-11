/*
Copyright 2014 Dominic Meiser

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
#ifndef SOURCE_SINK_HPP
#define SOURCE_SINK_HPP

#include <Amplitude.hpp>

namespace Lindblad {
namespace Detail {

/**
 * @brief Source for particles in a specific quantum level
 *
 * Particles are added to the system at a certain rate and particles are
 * removed from the system in any state at the same rate to ensure that
 * the overall number of particles is conserved (conservation of
 * trace[rho] == 1).
 * */
class SourceSink {
 public:
   /**
    * @brief Construct from rate and state.
    *
    * @param m     The quantum level in which atoms enter the system.
    * @param gamma The rate at which atoms enter the system.
    * */
  SourceSink(int m, double gamma) : m(m), gamma(gamma) {}
  void apply(int dim, const Amplitude *A, Amplitude *B) const;

 private:
  int m;
  double gamma;
};
}
}

#endif

