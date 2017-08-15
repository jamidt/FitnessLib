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

#include <string>
#include <sstream>
#include <vector>

namespace intrep {

template<class T>
std::string toString(const std::vector<T>& vec) {
    std::stringstream output;
    auto iter_vec = vec.cbegin();
    output << *iter_vec;
    for(++iter_vec; iter_vec != vec.cend(); ++iter_vec) {
        output << " " << *iter_vec;
    }
    return output.str();
}

template<std::size_t loci>
std::string str_vertex_int_bin(std::size_t vertex) {
    std::stringstream stream;
    stream << "(" << vertex << ", " << std::bitset<loci>(vertex).to_string() << ")";
    return stream.str();
}

template<std::size_t loci>
std::string str_vertex_bin(std::size_t vertex) {
    std::stringstream stream;
    stream << "(" << std::bitset<loci>(vertex).to_string() << ")";
    return stream.str();
}

} // end namespace intrep
