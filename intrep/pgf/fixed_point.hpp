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

#include <cstdlib>

namespace intrep {

template <class Container, class PGF>
class GenericFixpoint {
  private:
    PGF* _pgf;
    using calc_t = typename Container::calc_t;
    Container _c;

  public:
    GenericFixpoint(PGF* pgf) : _pgf(pgf) {}
    Container& calculateFixpoint() {
        _c.fill(0);
        for (std::size_t i = 0; i != 1000; ++i) {
            (*_pgf)(_c);
        }
        return &_c;
    }
};
}
