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
#include <gtest/gtest.h>
#include <GMRES.hpp>
#include <cmath>

void printVector(int dim, Amplitude* v) {
  for (int i = 0; i < dim; ++i) {
    std::cout << v[i] << std::endl;
  }
}

TEST(GMRES, Constructor) {
  GMRES *gmres = new GMRES(1); 
  ASSERT_TRUE(gmres != 0);
  delete gmres;
}

static void fDouble(int dim, const Amplitude *x, Amplitude *result,
                    void *ctx) {
  for (int i = 0; i < dim; ++i) {
    result[i] = 2.0 * x[i];
  }
}

TEST(GMRES, axpy) {
  GMRES gmres(2); 
  std::vector<Amplitude> y(2);
  y[0] = 2.0;
  y[1] = 3.0;
  std::vector<Amplitude> x(2);
  x[0] = -1.0;
  x[1] = 2.7;
  std::vector<Amplitude> result(2);
  double alpha = -2.4;
  gmres.axpy(alpha, &fDouble, &x[0], &y[0], &result[0], 0);
  EXPECT_FLOAT_EQ((alpha * 2.0 * x[0] + y[0]).real(), result[0].real());
}

TEST(GMRES, norm) {
  GMRES gmres(2); 
  std::vector<Amplitude> x(2, 1.3);
  double nrm = gmres.norm(2, &x[0]);
  EXPECT_FLOAT_EQ(sqrt(2.0 * 1.3 * 1.3), nrm);
}

TEST(GMRES, dot) {
  GMRES gmres(2); 
  std::vector<Amplitude> x(2, 1.3);
  double nrm = gmres.norm(2, &x[0]);
  Amplitude dot = gmres.dot(x.size(), &x[0], &x[0]);
  EXPECT_FLOAT_EQ(nrm * nrm, dot.real());
  EXPECT_FLOAT_EQ(0.0, dot.imag());
}

TEST(GMRES, solve) {
  GMRES gmres(2); 
  std::vector<Amplitude> b(2, 1.3);
//  b[0] = Amplitude(3.9, 2.7);
//  b[1] = Amplitude(1.4, 3.7);
  std::vector<Amplitude> x(2);
  gmres.solve(&fDouble, &b[0], &x[0], 0);
  EXPECT_FLOAT_EQ(0.5 * b[0].real(), x[0].real());
  EXPECT_FLOAT_EQ(0.5 * b[0].imag(), x[0].imag());
  EXPECT_FLOAT_EQ(0.5 * b[1].real(), x[1].real());
  EXPECT_FLOAT_EQ(0.5 * b[1].imag(), x[1].imag());
}
