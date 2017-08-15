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

#include "path_wright_fisher.hpp"
#include "wright_fisher.hpp"

namespace intrep {

template <std::size_t L, class IntType = long long, class RealType = double,
          template <std::size_t, class> class Container = CubeValues>
class Entropy {
  private:
    LinesOfDescentHash<L, IntType, RealType, Container> _lod;

  public:
    template <class S>
    Entropy(const Container<L, RealType>& fitness, RealType mutationrate,
            IntType popsize, IntType runs, IntType max_gen, S seed)
        : _lod(fitness, mutationrate, popsize, runs, max_gen, seed) {}
    Entropy(const Container<L, RealType>& fitness, RealType mutationrate,
            IntType popsize, IntType runs, IntType max_gen)
        : Entropy(fitness, mutationrate, popsize, runs, max_gen,
                  std::random_device{}()) {}
};

}  // end namespace intrep
