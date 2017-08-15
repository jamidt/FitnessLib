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

#include "generic_landscape.hpp"
#include <random>

/**
 * Old header. Do NOT use!!!!
 */

namespace intrep {

template <std::size_t L, template <std::size_t, class> class Container,
          class RandomDist>
class GenericHouseOfCards : public GenericLandscape<L, double, Container> {
  protected:
    std::mt19937_64 _rand_engine;
    RandomDist _dist;
    bool _escape;
    double _escape_value;

  public:
    GenericHouseOfCards()
        : _rand_engine(std::random_device{}()), _dist(), _escape(false) {}
    GenericHouseOfCards(typename RandomDist::param_type p)
        : GenericHouseOfCards(), _dist(p), _escape(false) {}
    GenericHouseOfCards(typename RandomDist::param_type p, double esc_value)
        : GenericHouseOfCards(p), _escape(true), _escape_value(esc_value) {}

    void createLandscape() {
        for (auto& f : this->_fitnes) {
            this->_fitnes = _dist(_rand_engine);
        }
        if (_escape) {
            this->_fitnes.back() = _escape_value;
        }
    }
    void createLandscape(typename RandomDist::param_type p) {
        _dist.param(p);
        _escape = false;
        createLandscape();
    }
    void createLandscape(typename RandomDist::param_type p,
                         double escape_value) {
        _escape_value = escape_value;
        _escape = true;
        createLandscape(p);
    }
};

template <std::size_t L,
          template <std::size_t, class> class Container = CubeValues>
using UniHouseOfCards =
    GenericHouseOfCards<L, Container, std::uniform_real_distribution<>>;

template <std::size_t L,
          template <std::size_t, class> class Container = CubeValues>
using ExpHouseOfCards =
    GenericHouseOfCards<L, Container, std::exponential_distribution<double>>;

/**@ingroup landscape
 * @brief House of cards landscape
 *
 * House of Cards model with uniform distribution
 *
 */
template <std::size_t L,
          template <std::size_t, class> class Container = CubeValues>
class UniformHouseOfCards : public GenericLandscape<L, double, Container> {
  private:
    std::mt19937_64 _rand_engine;
    std::uniform_real_distribution<> _randdist;
    double _escape;
    bool _esc;

  public:
    UniformHouseOfCards() : UniformHouseOfCards(0, 1) {}
    UniformHouseOfCards(double a, double b)
        : _rand_engine(std::random_device{}()) {
        createLandscape(a, b);
    }
    UniformHouseOfCards(double a, double b, double escape)
        : _rand_engine(std::random_device{}()), _escape(escape) {
        createLandscape(a, b, escape);
    }

    void createLandscape(double a, double b) {
        _esc = false;
        typename std::uniform_real_distribution<>::param_type new_param(a, b);
        _randdist.param(new_param);
        for (std::size_t i = 0; i != this->w.size(); ++i) {
            this->w[i] = _randdist(_rand_engine);
        }
    }
    void createLandscape(double a, double b, double escape) {
        createLandscape(a, b);
        _esc = true;
        this->w.back() = escape;
    }

    void resetLandscape() {
        for (auto& x : this->w) {
            x = _randdist(_rand_engine);
        }
        if (_esc) {
            this->w.back() = _escape;
        }
    }
};
}
