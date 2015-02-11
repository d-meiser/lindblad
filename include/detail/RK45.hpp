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
#ifndef RK45_HPP
#define RK45_HPP

#include <vector>

#include <Integrator.hpp>

namespace Lindblad {
namespace Detail {

class RK45 : public Integrator {
 public:
  RK45(int dim, double time, const double* state,
      void (*f)(double* x, double* y, double t, void* ctx));
  ~RK45();

 private:
  std::vector<double> y;
  std::vector<std::vector<double> > ks;
  std::vector<double> work;
  std::vector<double> sol4;
  std::vector<double> sol5;

  static const double cs[6];
  static const double as[25];
  static const double b4[6];
  static const double b5[6];

  void computeWork(int row, double dt);
  double errorEstimate(const std::vector<double>& x, const std::vector<double>& y) const;
  virtual const double* getCurrentState() const;
  virtual void advance(double* t, double* dt, void* ctx);
  virtual RK45* makeCopy() const;
};
}
}

#endif
