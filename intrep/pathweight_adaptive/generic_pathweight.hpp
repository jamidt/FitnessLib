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
#include "../pathways.hpp"

#include <unordered_map>
#include <utility>
// NOTE: Boost has to be installed. No std::hash for std::array :`(
#include "boost/functional/hash.hpp"

#include <bitset>
#include <sstream>
#include <string>

namespace intrep {

// Container for the two verticies with pair interaction (i.e. possible
// steps).
using pair_inter = std::pair<std::size_t, std::size_t>;

/*
 * Print pair of loci in the form
 * (int(i), binary(i)) -> (int(j), binary(j))
 */
template <std::size_t loci>
std::string str_pair_bin(pair_inter p) {
    std::stringstream stream;
    stream << intrep::str_vertex_bin<loci>(p.first) << " -> "
           << intrep::str_vertex_bin<loci>(p.second);
    return stream.str();
}

/*
 * Print pair of loci in the form
 * int(i) -> int(j)
 */
template <std::size_t loci>
std::string str_pair_int(pair_inter p) {
    std::stringstream stream;
    stream << p.first << " -> " << p.second;
    return stream.str();
}

/*
 * Class saves the selection coefficent, fixation probability and step
 * probability at a vertex.
 */
template <class T>
struct values {
    T fixation_prob;
    T sel_coeff;
    T step_prob;
    values()
        : fixation_prob(0),
          sel_coeff(0),
          step_prob(0) {}  // Default constr for unordered_map
    values(T sel, T fix, T step)
        : fixation_prob(fix), sel_coeff(sel), step_prob(step) {}
};

/**@ingroup adaptivepathweight
 * @brief Generic class for accessible pathweights
 *
 * TODO: This is suboptimal. The precision should not be a template
 * parameter of the class. In fact, not even the number of loci should
 * be a template parameter here, because it is not needed here.
 *
 * Generic class for pathweight calculation of adaptive walks on a
 * single peaked landscape. See Weinreich 2006.
 * @tparam loci
 * @tparam T
 */
template <std::size_t loci, class T = double>
class GenericAccessiblePathweight {
  protected:
    /// Container for the weights of each step
    std::unordered_map<pair_inter, values<T>, boost::hash<pair_inter>>
        _pair_value;

  public:
    // Get values for a pair
    values<T> get(pair_inter p) { return _pair_value.find(p); }

    /**@brief Get pathprobability of a given path
     *
     * @tparam P Container storing the pathway. This datastructure must
     * implement a const_iterator.
     * @param[in] path Pathway for which the weight should be
     * calculated. Note that the path is given from the starting point
     * to the final (end) type.
     */
    template <class P>
    T getPathProb(const P& path) {
        auto iter1 = path.cbegin();
        auto iter2 = path.cbegin();
        ++iter2;
        T path_prob = static_cast<T>(1);
        for (; iter2 != path.cend(); ++iter2, ++iter1) {
            path_prob =
                path_prob *
                _pair_value.at(std::make_pair(*iter1, *iter2)).step_prob;
        }
        return path_prob;
    }

    using iterator =
        typename std::unordered_map<pair_inter, values<T>,
                                    boost::hash<pair_inter>>::iterator;
    using const_iterator =
        typename std::unordered_map<pair_inter, values<T>,
                                    boost::hash<pair_inter>>::const_iterator;

    ///@name Iterator
    ///@{
    iterator begin() { return _pair_value.begin(); }
    iterator end() { return _pair_value.end(); }
    const_iterator cbegin() { return _pair_value.cbegin(); }
    const_iterator cend() { return _pair_value.cend(); }
    ///@}
};
}  // end namespace intrep
