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
#include <RK45.hpp>

RK45::RK45(int dim, double time, const double* state,
           void (*f)(double* x, double* y, double t, void* ctx))
    : Integrator(dim, time, state, f) {}
const double* RK45::getCurrentState() const {}
void RK45::advance(double* t, double* dt, void* ctx) {}
RK45* RK45::makeCopy() const {}
