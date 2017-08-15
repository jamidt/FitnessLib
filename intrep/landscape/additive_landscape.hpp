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

/**
 * Additive landscape. The Constructor and the function
 * createLandscape take the values of the wilttype, all single
 * mutation types and the escaping type as argument.
 */

#include "../math.hpp"
#include "generic_landscape.hpp"
#include <utility>

namespace intrep {
/**@ingroup landscape
 * @brief Generic additive landscape
 *
 * GenericAdditiveLandscape is constructed from the fitness values of the
 * first L single mutants.
 */
template <std::size_t loci, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class GenericAdditiveLandscape : public GenericLandscape<loci, T, Container> {
  private:
    template <class U>
    void setStartHelper(U&& set_type) {
        this->_fitness[1] = set_type;
    }
    template <class Head, class... Tail>
    void setStartHelper(Head&& set_type, Tail&&... tail) {
        size_t l = 1 << sizeof...(tail);
        this->_fitness[l] = set_type;
        setStartHelper(std::forward<Tail>(tail)...);
    }

  public:
    GenericAdditiveLandscape() = default;
    /// Constructor
    // w(10..0), w(010..0), .., w(0..01)
    template <class... Args>
    GenericAdditiveLandscape(Args&&... args) {
        createLandscape(std::forward<Args>(args)...);
    }
    // Same syntax as with the constructor
    template <class... Arg>
    void createLandscape(Arg&&... args) {
        if (sizeof...(args) == loci) {
            this->_fitness.front() = 0;
            setStartHelper(0, std::forward<Arg>(args)...);
            for (size_t vertex = 0; vertex != verticies(loci); ++vertex) {
                if (hammingweight(vertex) > 1) {
                    for (size_t n = 0; n != loci; ++n) {
                        size_t desc = 1 << n;
                        if ((vertex & desc) != 0) {
                            this->_fitness[vertex] += this->_fitness[desc];
                        } else {
                            continue;
                        }
                    }
                }
            }
        } else {
            throw std::length_error(
                "The L (=loci) first fitness values must be given\n");
        }
    }
};

/**@ingroup landscape
 * @brief Additive landscape with offset.
 *
 * Landscape is given by a simple
 * additive landscape \f$F'(\sigma)\f$ and a offset \f$c\f$.
 *
 * \f$F(\sigma)=F'(\sigma)+c\f$
 */
template <std::size_t loci, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class AdditiveOffsetLandscape
    : public GenericAdditiveLandscape<loci, T, Container> {
  public:
    template <class Offset, class... Tail>
    void createLandscape(Offset&& offset,
                         Tail&&... arg) {  // Why no r-value for T?!??!!???
        if (sizeof...(arg) == loci) {
            GenericAdditiveLandscape<loci, T, Container>::createLandscape(
                std::forward<Tail>(arg)...);
            for (auto& x : this->_fitness) {
                x += offset;
            }
        } else {
            std::stringstream err_string;
            err_string << loci + 1
                       << " values must be given in the following order: "
                          "w(0..0), w(10..0), ..., w(0..01). Arguments given "
                       << sizeof...(arg) + 1 << "\n";
            throw std::length_error(err_string.str());
        }
    }
    /**@brief Constructor
     *
     * First argument is the offset value. All other values provide the
     * fitness values in the order \f$(10\dots 0), (010\dots 0), \dots,
     * (0\dots 01)\f$.
     */
    template <class Offset, class... Tail>
    AdditiveOffsetLandscape(Offset&& offset, Tail&&... arg) {
        createLandscape(std::forward<Offset>(offset),
                        std::forward<Tail>(arg)...);
    }

    AdditiveOffsetLandscape() = default;
};

/**@ingroup landscape
 * @brief Additive subcritical landscape
 *
 * Additive fitness landscape, begining with an offset value at (00...0)
 * and all fitness values subcritical but the L mutant (11...1) which is
 * given by a supercritical fitness value.
 */
template <std::size_t loci, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class AdditiveSubcriticalEscapeLandscape
    : public AdditiveOffsetLandscape<loci, T, Container> {
  public:
    template <class Escape, class Offset, class... Tail>
    void createLandscape(Escape&& escape, Offset&& offset, Tail&&... arg) {
        AdditiveOffsetLandscape<loci, T, Container>::createLandscape(
            std::forward<Offset>(offset), std::forward<Tail>(arg)...);
        this->_fitness.back() = escape + this->_fitness.front();
    }
    template <class Escape, class Offset, class... Tail>
    AdditiveSubcriticalEscapeLandscape(Escape&& escape, Offset&& offset,
                                       Tail&&... arg) {
        if (sizeof...(arg) == loci) {
            createLandscape(std::forward<Escape>(escape),
                            std::forward<Offset>(offset),
                            std::forward<Tail>(arg)...);
        } else {
            std::stringstream err_string;
            err_string << loci + 1 << " values must be given\n";
            throw std::length_error(err_string.str());
        }
    }
};

/**@ingroup landscape
 * @brief Geometric additive landscape
 *
 * Landscape is constructed
 * \f[
 * F(\sigma)=\sum_{k=1}^{L}f_k\sigma_\k
 * \f]
 * with \f$f_k = ab^k \f$
 */
template <std::size_t loci, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class GeometricAdditiveLandscape
    : public GenericAdditiveLandscape<loci, T, Container> {
  public:
    /**@brief Constructor
     */
    GeometricAdditiveLandscape(T a, T b) { createLandscape(a, b); }
    GeometricAdditiveLandscape() = default;
    void createLandscape(T a, T b) {
        createLandscape_imp(a, b, std::make_index_sequence<loci>{});
    }

  private:
    template <std::size_t... I>
    void createLandscape_imp(T a, T b, std::index_sequence<I...>) {
        GenericAdditiveLandscape<loci, T, Container>::createLandscape(
            makeFitness_impl(a, b, I)...);
    }
    T makeFitness_impl(T a, T b, std::size_t k) {
        return a * intrep::pow(b, static_cast<long>(k + 1));
    }
};

/**@ingroup landscape
 * @brief Geometric additive landscape with offset
 *
 * Landscape is constructed
 * \f[
 * F(\sigma)=\sum_{k=1}^{L}f_k\sigma_\k + f_0
 * \f]
 * with \f$f_k = ab^k \f$
 */
template <std::size_t loci, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class GeometricAdditiveOffsetLandscape
    : public GeometricAdditiveLandscape<loci, T, Container> {
  public:
    GeometricAdditiveOffsetLandscape(T a, T b, T offset) {
        createLandscape(a, b, offset);
    }
    GeometricAdditiveOffsetLandscape() = default;
    void createLandscape(T a, T b, T offset) {
        GeometricAdditiveLandscape<loci, T, Container>::createLandscape(a, b);
        for (auto& f : this->_fitness) {
            f = f + offset;
        }
    }
};

/**@ingroup landscape
 * @brief Geometric offset landscape with set escape type
 */
template <std::size_t loci, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class GeometricAdditiveOffsetEscapeLandscape
    : public GeometricAdditiveOffsetLandscape<loci, T, Container> {
  public:
    GeometricAdditiveOffsetEscapeLandscape(T a, T b, T offset, T escape) {
        createLandscape(a, b, offset, escape);
    }
    GeometricAdditiveOffsetEscapeLandscape() = default;
    void createLandscape(T a, T b, T offset, T escape) {
        GeometricAdditiveOffsetLandscape<loci, T, Container>::createLandscape(
            a, b, offset);
        this->_fitness.back() = escape;
    }
};

}  // end namespace intrep
