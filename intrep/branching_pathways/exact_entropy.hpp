/**
 * Copyright © 2017 Jan Schmidt
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../landscape.hpp"
#include "../math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>

#include "../../MP/mpreal.h"
#include "../pathways.hpp"
#include "../pgf/hypercube_pgf.hpp"
#include "path.hpp"

namespace intrep {
template <std::size_t L, class T, class F = T,
          template <std::size_t, class> class Container = CubeValues>
class ExactBranchingDistSingle {
  private:
    const GenericLandscape<L, F>& _fitness;
    F _mutationrate;
    Pgf01<L, T, F, Container, Container> _pgf;
    Pathways01<L> _pathways;
    using ProbDist_t = std::array<T, factorial(L)>;
    ProbDist_t _path_probabilities;
    const int _generations;

  public:
    ExactBranchingDistSingle(const GenericLandscape<L, F>& fitness,
                             F mutationrate, int generations)
        : _fitness(fitness),
          _mutationrate(mutationrate),
          _pgf(_fitness, mutationrate),
          _pathways(),
          _generations(generations) {}

    const ProbDist_t& getProbDist() const { return _path_probabilities; }
    const Pathways01<L>& getPathways() const { return _pathways; }

    void createDist() {
        for (std::size_t n = 0; n != _pathways.size(); ++n) {
            BranchingTimeDist01<L, T, Container> path_d;
            for (int gen = 0; gen != _generations; ++gen) {
                path_d.nextTimestep(_pgf, _pathways[n]);
            }
            _path_probabilities.at(n) = path_d.getCProb();
        }
        normalize(_path_probabilities);
    }

    const ProbDist_t& getPathDist() const { return _path_probabilities; }
    using iterator = typename ProbDist_t::iterator;
    using const_iterator = typename ProbDist_t::const_iterator;
    iterator begin() { return _path_probabilities.begin(); }
    const_iterator begin() const { return _path_probabilities.cbegin(); }
    const_iterator cbegin() const { return _path_probabilities.cbegin(); }
    iterator end() { return _path_probabilities.end(); }
    const_iterator end() const { return _path_probabilities.cend(); }
    const_iterator cend() const { return _path_probabilities.cend(); }

    T calcPathEntropy() const {
        T entropy = std::accumulate(
            _path_probabilities.cbegin(), _path_probabilities.cend(),
            static_cast<T>(0),
            [](T& ent, const T& p) { return ent - p * intrep::log(p); });
        return entropy;
    }
    T calcKLDivergence() const {
        T uniform = 1 / static_cast<T>(_fitness.size());
        T kldiv = std::accumulate(_path_probabilities.cbegin(),
                                  _path_probabilities.cend(), static_cast<T>(0),
                                  [=](T& ent, const T& p) {
                                      return ent + p * intrep::log(p / uniform);
                                  });
        return kldiv;
    }
};

template <std::size_t L,
          template <std::size_t, class> class Container = CubeValues>
using MprealExactBranchingDistSingle =
    ExactBranchingDistSingle<L, mpfr::mpreal, mpfr::mpreal, Container>;
}
