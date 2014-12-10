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
#ifndef GMRES_HPP
#define GMRES_HPP

#include <vector>
#include <Amplitude.hpp>

class GMRES {
public:
  GMRES(int dim);
  void solve(void (*Ax)(const Amplitude *x, Amplitude *y, void *ctx),
             const Amplitude *rhs, Amplitude *x);

private:
  std::vector<Amplitude> y;
  static const int m = 30;
};

#endif

