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
#include "generic_landscape.hpp"
#include <random>

namespace intrep {


/** House of Cards landscape
 *
 * The template parameter T is for the container and double must be
 * castable to T.
 */
template <std::size_t L, class T, template <std::size_t, class> class Container,
          class RandomDist>
class GenericHouseOfCards : public GenericLandscape<L, T, Container> {
  private:
    std::mt19937_64 _rand_engine;
    RandomDist _dist;
    bool _escape;
    T _escape_value;

  public:
    GenericHouseOfCards()
        : _rand_engine(std::random_device{}()), _dist(), _escape(false) {
        createLandscape();
    }
    GenericHouseOfCards(typename RandomDist::param_type p)
        : _rand_engine(std::random_device{}()), _dist(p), _escape(false) {
        createLandscape();
    }
    GenericHouseOfCards(typename RandomDist::param_type p, T esc_value)
        : GenericHouseOfCards(p) {
        _escape_value = esc_value;
        _escape = true;
        createLandscape();
    }

    void setParam(typename RandomDist::param_type p) {
        _escape = false;
        _dist.param(p);
    }
    void setParam(typename RandomDist::param_type p, T escape) {
        _escape = true;
        _escape_value = escape;
        _dist.param(p);
    }

    void createLandscape() {
        for (auto& f : this->_fitness) {
            f = static_cast<T>(_dist(_rand_engine));
        }
        if (_escape) {
            this->_fitness.back() = _escape_value;
        }
    }
    void createLandscape(typename RandomDist::param_type p) {
        _dist.param(p);
        _escape = false;
        for(auto& f: this->_fitness) {
            f = static_cast<T>(_dist(_rand_engine));
        }
    }
    void createLandscape(typename RandomDist::param_type p,
                         double escape_value) {
        _escape_value = escape_value;
        _escape = true;
        _dist.param(p);
        // Unfortunately, the compiler does something stupid when using
        // optimizations. The value true for _escape is not set before
        // calculating the random parts. Hence, this has to be done here
        // manually... :`(
        for(auto& f: this->_fitness) {
            f = static_cast<T>(_dist(_rand_engine));
        }
        this->_fitness.back() = _escape_value;
    }
};

template <std::size_t L, class T = double,
          template <std::size_t, class> class Container = CubeValues>
using UniHouseOfCards =
    GenericHouseOfCards<L, T, Container, std::uniform_real_distribution<>>;

template <std::size_t L, class T = double,
          template <std::size_t, class> class Container = CubeValues>
using ExpHouseOfCards =
    GenericHouseOfCards<L, T, Container, std::exponential_distribution<double>>;

// Rough Mt Fuji landscape
template <std::size_t L, template <std::size_t, class> class Container,
          class RandomDist>
class GenericRMF : public GenericLandscape<L, double, Container> {
  private:
    std::mt19937_64 _rand_engine;
    RandomDist _dist;
    double _slope;

  public:
    GenericRMF(double slope)
        : _rand_engine(std::random_device{}()), _dist(), _slope(slope) {
        createLandscape();
    }
    GenericRMF(typename RandomDist::param_type p, double slope = static_cast<double>(1))
        : _rand_engine(std::random_device{}()), _dist(p), _slope(slope) {
        createLandscape();
    }
    GenericRMF() : GenericRMF(static_cast<double>(1)) { createLandscape(); }

    void setParam(typename RandomDist::param_type p) {
        _dist.param(p);
    }
    void setParam(typename RandomDist::param_type p, double slope) {
        _slope = slope;
        _dist.param(p);
    }

    void createLandscape() {
        for (std::size_t vertex = 0; vertex != this->_fitness.size(); ++vertex) {
            this->_fitness[vertex] =
                _slope * hammingweight(vertex) + _dist(_rand_engine);
        }
    }
    void createLandscape(typename RandomDist::param_type p) {
        _dist.param(p);
        for (std::size_t vertex = 0; vertex != this->_fitness.size(); ++vertex) {
            this->_fitness[vertex] =
                _slope * hammingweight(vertex) + _dist(_rand_engine);
        }
    }
    void createLandscape(typename RandomDist::param_type p, double slope) {
        _slope = slope;
        _dist.param(p);
        for (std::size_t vertex = 0; vertex != this->_fitness.size(); ++vertex) {
            this->_fitness[vertex] =
                _slope * hammingweight(vertex) + _dist(_rand_engine);
        }
    }
};

template <std::size_t L,
          template <std::size_t, class> class Container = CubeValues>
using UniformRMF = GenericRMF<L, Container, std::uniform_real_distribution<>>;

template <std::size_t L,
          template <std::size_t, class> class Container = CubeValues>
using ExpRMF = GenericRMF<L, Container, std::exponential_distribution<>>;

}  // end namespace intrep

template<std::size_t L, template<std::size_t, class> class Container, class RandomDist>
std::ostream& operator<<(std::ostream& out, const intrep::GenericRMF<L, Container, RandomDist>& fitness) {
    for(std::size_t vertex = 0; vertex != fitness.size(); ++vertex) {
        out << vertex << ' ' << fitness[vertex] << '\n';
    }
    return out;
}
