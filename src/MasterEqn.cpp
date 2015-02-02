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
#include <vector>

struct MasterEqn::Impl {
  Impl(int d) : dim(d) {}
  ~Impl() {}
  void apply(const Amplitude* A, Amplitude* B) const {
    std::fill(B, B + dim * dim, 0);
    for (std::vector<Coupling>::const_iterator c = couplings.begin();
         c != couplings.end(); ++c) {
      c->apply(dim, A, B);
    }
    for (std::vector<Decay>::const_iterator d = decays.begin();
         d != decays.end(); ++d) {
      d->apply(dim, A, B);
    }
    for (std::vector<SourceSink>::const_iterator s = sourceSinks.begin();
         s != sourceSinks.end(); ++s) {
      s->apply(dim, A, B);
    }
  }
  int getDim() const { return dim; }

  int dim;
  std::vector<Coupling> couplings;
  std::vector<Decay> decays;
  std::vector<SourceSink> sourceSinks;
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
  impl->couplings.push_back(Coupling(m, n, a));
}

void MasterEqn::addDecay(int into, int outOf, double gamma) {
  impl->decays.push_back(Decay(into, outOf, gamma));
}

void MasterEqn::addSourceSink(int into, double gamma) {
  impl->sourceSinks.push_back(SourceSink(into, gamma));
}

void MasterEqn::apply(const Amplitude* A, Amplitude *B) const {
  impl->apply(A, B);
}

int MasterEqn::getDim() const {
  return impl->getDim();
}

