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
#include <iostream>
#include <vector>

#define STREAMKEYVALUE(stream, key) \
    do { (stream) << "# " << #key << " = " << (key) << '\n'; } while(false);

namespace intrep {
namespace io {
enum class loci : std::size_t {
    integer = 0x1,
    binary = 0x2,
    integer_binary = 0x4,
};

loci operator|(loci x, loci y) {
    return static_cast<loci>(static_cast<std::size_t>(x) |
                             static_cast<std::size_t>(y));
}

std::ostream& operator<<(std::ostream& out, const loci& l) {
    out << static_cast<std::size_t>(l);
    return out;
}

}
}
// Printing lines of descent (that is, the loci in integer
// representation)
template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& l) {
    auto iter = l.cbegin();
    auto end = l.cend();
    --end;
    for (; iter != end; ++iter) {
        out << *iter << ' ';
    }
    out << *end;
    return out;
}
