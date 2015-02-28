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
#include <MasterEqn.hpp>
#include <Coupling.hpp>
#include <Decay.hpp>
#include <SourceSink.hpp>
#include <GeneralDecayOperator.hpp>
#include <vector>

namespace Lindblad {

struct MasterEqn::Impl {
  Impl(int d) : dim(d) {}
  ~Impl() {}
  void apply(const Amplitude* A, Amplitude* B) const {
    std::fill(B, B + dim * dim, 0);
    for (std::vector<Detail::Coupling>::const_iterator c = couplings.begin();
         c != couplings.end(); ++c) {
      c->apply(dim, A, B);
    }
    for (std::vector<Detail::Decay>::const_iterator d = decays.begin();
         d != decays.end(); ++d) {
      d->apply(dim, A, B);
    }
    for (std::vector<Detail::SourceSink>::const_iterator s =
             sourceSinks.begin();
         s != sourceSinks.end(); ++s) {
      s->apply(dim, A, B);
    }
    for (std::vector<Detail::GeneralDecayOperator>::const_iterator d =
             generalDecayOperators.begin();
         d != generalDecayOperators.end(); ++d) {
      d->apply(dim, A, B);
    }
  }
  int getDim() const { return dim; }
  void getEnergyLevels(Amplitude* omegas) const {
    std::fill(omegas, omegas + dim, 0);
    for (std::vector<Detail::Coupling>::const_iterator c = couplings.begin();
         c != couplings.end(); ++c) {
      if (c->isDiagonal()) {
        omegas[c->getRow()] += c->getCouplingStrength();
      }
    }
  }

  int dim;
  std::vector<Detail::Coupling> couplings;
  std::vector<Detail::Decay> decays;
  std::vector<Detail::SourceSink> sourceSinks;
  std::vector<Detail::GeneralDecayOperator> generalDecayOperators;
};

MasterEqn::MasterEqn(int d) { impl = new Impl(d); }

MasterEqn::MasterEqn(const MasterEqn& other)
    : impl(new MasterEqn::Impl(*other.impl)) {}

MasterEqn& MasterEqn::operator=(const MasterEqn& rhs) {
  if (this != &rhs) {
    delete impl;
    impl = new MasterEqn::Impl(*rhs.impl);
  }
  return *this;
}

MasterEqn::~MasterEqn() {
  delete impl;
}

void MasterEqn::addCoupling(int m, int n, Amplitude a) {
  impl->couplings.push_back(Detail::Coupling(m, n, a));
}

void MasterEqn::addDecay(int into, int outOf, double gamma) {
  impl->decays.push_back(Detail::Decay(into, outOf, gamma));
}

void MasterEqn::addSourceSink(int into, double gamma) {
  impl->sourceSinks.push_back(Detail::SourceSink(into, gamma));
}

void MasterEqn::addGeneralDecayOperator(SparseMatrix lambda) {
  impl->generalDecayOperators.push_back(Detail::GeneralDecayOperator(lambda));
}

void MasterEqn::apply(const Amplitude* A, Amplitude *B) const {
  impl->apply(A, B);
}

int MasterEqn::getDim() const {
  return impl->getDim();
}

void MasterEqn::getEnergyLevels(Amplitude *omegas) const {
  impl->getEnergyLevels(omegas);
}

static void buildIdentityMatrix(int d, Amplitude* m) {
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
      m[i * d + j] = (i == j) ? Amplitude(1) : Amplitude(0);
    }
  }
}

static void transposeMatrix(int d, Amplitude* m) {
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
      Amplitude tmp = m[i * d + j];
      m[i * d + j] = m[j * d + i];
      m[j * d + i] = tmp;
    }
  }
}

void MasterEqn::buildMatrix(Amplitude* matrix) const {
  buildTransposedMatrix(matrix);
  int dim = getDim();
  transposeMatrix(dim * dim, &matrix[0]);
}

void MasterEqn::buildTransposedMatrix(Amplitude* matrix) const {
  int dim = getDim();
  std::vector<Amplitude> identityMatrix(dim * dim * dim * dim);
  buildIdentityMatrix(dim * dim, &identityMatrix[0]);
  for (int i = 0; i < dim * dim; ++i) {
    apply(&identityMatrix[i * dim * dim], &matrix[i * dim * dim]);
  }
}

}
