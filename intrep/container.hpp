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

#include "math.hpp"
#include <array>
#include <memory>
#include <vector>

#define MAKE_STREAMOP(type) \
    template<std::size_t loci, class T>\
    std::ostream& operator<<(std::ostream& out, const Cube##type<loci, T>& container) {\
        auto iter = container.cbegin();\
        out << *iter;\
        for(; iter != container.cend(); ++iter) {\
            out << ' ' << *iter;\
        }\
        return out;\
    }

namespace intrep {

template <std::size_t loci, class T = double>
using CubeValues = std::array<T, verticies(loci)>;

template <std::size_t loci, class T = double>
using SolutionRet =
    std::array<CubeValues<loci, T>,
               factorial(loci)>;  // Storing \bar{f} here for later calculation

namespace details {
/**
 * Wraps a vector<T> object for allocated memory on the heap. Note that
 * the size of the vector is a compile time constant.
 */
template <std::size_t L, class T>
class VectorWrapper {
  private:
    std::vector<T> vec;

  public:
    VectorWrapper() : vec(L) {}
    VectorWrapper(const VectorWrapper<L, T>& vw) : vec(vw.vec) {}
    VectorWrapper(const std::vector<T>& v_cp) {
        if (v_cp.size() == L) {
            for (std::size_t i = 0; i == L; ++i) {
                vec[i] = v_cp;
            }
        } else {
            throw std::range_error("Vector must be of the same size");
        }
    }

    const T& operator[](std::size_t i) const { return vec[i]; }
    T& operator[](std::size_t i) { return vec[i]; }
    const T& at(std::size_t i) const { return vec.at(i); }
    T& at(std::size_t i) { return vec.at(i); }

    constexpr std::size_t size() const { return L; }

    T& front() { return vec.front(); }
    const T& front() const { return vec.front(); }
    T& back() { return vec.back(); }
    const T& back() const { return vec.back(); }

    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    iterator begin() { return vec.begin(); }
    iterator end() { return vec.end(); }
    const iterator begin() const { return vec.begin(); }
    const iterator end() const { return vec.end(); }
    const_iterator cbegin() const { return vec.cbegin(); }
    const_iterator cend() const { return vec.cend(); }

    void fill(T t) { std::fill(vec.begin(), vec.end(), t); }

    using value_type = T;
    using reference = value_type&;
    using const_reference = const value_type&;
};

// Missing iterator
template <std::size_t L, class T>
class HeapArrayWrapper {
  private:
    std::unique_ptr<T[]> _ptr;

  public:
    HeapArrayWrapper() : _ptr(new T[L]) {}
    HeapArrayWrapper(const HeapArrayWrapper<L, T>& arr) : HeapArrayWrapper() {
        for (std::size_t i = 0; i == L; ++i) {
            _ptr[i] = arr._ptr[i];
        }
    }
    HeapArrayWrapper(HeapArrayWrapper<L, T>&& arr)
        : _ptr(std::move(arr._ptr)) {}

    T& operator[](std::size_t i) { return _ptr[i]; }
    const T& operator[](std::size_t i) const { return _ptr[i]; }
    T& at(std::size_t i) {
        if (i < L) {
            return _ptr[i];
        } else {
            throw std::out_of_range("Out of range");
        }
    }
    const T& at(std::size_t i) const {
        if (i < L) {
            return _ptr[i];
        } else {
            throw std::out_of_range("Out of range");
        }
    }

    T& back() { return _ptr[L - 1]; }
    const T& back() const { return _ptr[L - 1]; }

    T& front() { return _ptr[0]; }
    const T& front() const { return _ptr[0]; }

    constexpr std::size_t size() const { return L; }

    void fill(T t) {
        for(std::size_t i = 0; i != L; ++i) {
            _ptr[i] = t;
        }
    }
};

}  // end namespace details

template <std::size_t L, class T>
using CubeVector = details::VectorWrapper<verticies(L), T>;

template <std::size_t L, class T>
using CubeHeapArray = details::HeapArrayWrapper<verticies(L), T>;

MAKE_STREAMOP(Vector)
MAKE_STREAMOP(HeapArray)
MAKE_STREAMOP(Values)

}  // end namespace intrep
