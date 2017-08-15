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

#include "../MP/mpreal.h"
#include <algorithm>
#include <cmath>
#include <type_traits>

/**@defgroup math Mathfunctions
 * @brief Here some compiletime functions and selections for types
 * and mpfr::mpreal values are done
 */

namespace intrep {
/**
 * Normalize container
 */
template <class Container>
void normalize(Container& container) {
    typename Container::value_type normalization =
        std::accumulate(container.cbegin(), container.cend(),
                        static_cast<typename Container::value_type>(0));
    std::for_each(
        container.begin(), container.end(),
        [=](typename Container::value_type& x) { x = x / normalization; });
}

/**@ingroup math
 * @brief Caculates the number of verticies
 *
 * Function calculates the number of verticies for given loci during compile
 * time (i.e. \f$2^{loci}\f$)
 */
constexpr std::size_t verticies(std::size_t loci) {
    return loci > 0 ? 2 * verticies(loci - 1) : 1;
}

/**@ingroup math
 * @brief \f$n!\f$
 *
 * This calculates the factorial during compile time
 */
constexpr std::size_t factorial(std::size_t n) {
    return n > 0 ? n * factorial(n - 1) : 1;
}

// Calculate the Hamming weight
// TODO: use __buildin_popcountl() ?
// http://wm.ite.pl/articles/sse-popcount.html
// https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSet64
// https://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
template <class T>
struct select_int_t {};

template <>
struct select_int_t<std::uint64_t> {
    static uint64_t popcount(uint64_t vertex) {
        const uint64_t m1 = UINT64_C(0x5555555555555555);
        const uint64_t m2 = UINT64_C(0x3333333333333333);
        const uint64_t m4 = UINT64_C(0x0F0F0F0F0F0F0F0F);
        const uint64_t h01 = UINT64_C(0x0101010101010101);

        vertex -= (vertex >> 1) & m1;
        vertex = (vertex & m2) + ((vertex >> 2) & m2);
        vertex = (vertex + (vertex >> 4)) & m4;
        return (vertex * h01) >> 56;
    }
};

template <class T>
T hammingweight(T vertex) {
    return select_int_t<T>::popcount(vertex);
}

template <class T>
struct func {
    static T exp(const T& x) { return std::exp(x); }
    template <class E>
    static T pow(const T& x, const E& e) {
        return std::pow(x, e);
    }
    static T abs(const T& x) { return std::abs(x); }
    static T log(const T& x) { return std::log(x); }
};

template <>
struct func<mpfr::mpreal> {
    static mpfr::mpreal exp(const mpfr::mpreal& x) { return mpfr::exp(x); }
    template <class E>
    static mpfr::mpreal pow(mpfr::mpreal x, E e) {
        return mpfr::pow(x, e);
    }
    static mpfr::mpreal abs(const mpfr::mpreal& x) { return mpfr::abs(x); }
    static mpfr::mpreal log(const mpfr::mpreal& x) { return mpfr::log(x); }
};

/**@ingroup math
 * @brief \f$e^x\f$
 *
 * This calculates the exponential for a given type. Note this is done
 * only for selecting std::exp for double and mpfr::exp for
 * mpfr::mpreal.
 */
template <class T>
T exp(const T& x) {
    return func<T>::exp(x);
}

/**@ingroup math
 * @brief \f$x^y\f$
 *
 * This calculates the exponential for a given type. Note this is done
 * only for selecting std::pow for types and mpfr::pow for
 * mpfr::mpreal.
 */
template <class T, class E>
T pow(T x, E y) {
    return func<T>::pow(x, y);
}

/**@ingroup math
 * @brief \f$|x|\f$
 *
 * This calculates the exponential for a given type. Note this is done
 * only for selecting std::abs for types and mpfr::abs for
 * mpfr::mpreal.
 */
template <class T>
T abs(const T& x) {
    return func<T>::abs(x);
}

template<class T>
T log(const T& x) {
    return func<T>::log(x);
}

}  // end namespace intrep
