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
#include <SourceSink.hpp>

void SourceSink::apply(int dim, const Amplitude *rhoIn,
                       Amplitude *rhoOut) const {
  for (int i = 0; i < dim * dim; ++i) {
    rhoOut[i] -= gamma * rhoIn[i];
  }
  Amplitude trace(0, 0);
  for (int i = 0; i < dim; ++i) {
    trace += rhoIn[i + i * dim];
  }
  rhoOut[m + m * dim] += trace * gamma;
}

