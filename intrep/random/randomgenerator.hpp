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

#include "../math.hpp"
#include <random>

namespace intrep {
class BaseRandomGenerator {
  protected:
    std::mt19937_64 _rand_eng;

  public:
    BaseRandomGenerator() : _rand_eng(std::random_device{}()) {}
};

template <class Derived>
class RMFGeneric : public BaseRandomGenerator {
  protected:
    double _slope;

  public:
    RMFGeneric(double slope) : BaseRandomGenerator(), _slope(slope) {}
    double operator()(std::size_t vertex) {
        return hammingweight(vertex) * _slope +
               static_cast<Derived*>(this)->dist(this->_rand_eng);
    }
};

class RMFUniform : public RMFGeneric<RMFUniform> {
  private:
    std::uniform_real_distribution<> dist;

  public:
    RMFUniform(double slope) : RMFUniform(slope, 0, 1) {}
    RMFUniform(double slope, double min, double max)
        : RMFGeneric(slope), dist(min, max) {}
};

class RMFExponential : public RMFGeneric<RMFExponential> {
  private:
    std::exponential_distribution<> _dist;

  public:
    RMFExponential(double slope, double lambda)
        : RMFGeneric(slope), _dist(lambda) {}
};

}  // end namespace intrep
