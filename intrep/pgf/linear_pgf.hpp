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

#include <array>
#include "../math.hpp"
#include <type_traits>
#include <functional>

#include <iostream>

namespace intrep {
    template<std::size_t L, class T = double>
    using container_linear = std::array<T, L>;

    template<std::size_t L, class T>
    class GenericLinearPGF {
        public:
            void operator()(container_linear<L, T>&) = 0;
    };

    template<std::size_t L, class T = double>
    class BinarySplittingPGF : GenericLinearPGF<L, T> {
        private:
            container_linear<L, T> _tmp;
            const container_linear<L, T>& _birth;
            T _mutation;
        public:
            BinarySplittingPGF(const container_linear<L, T>& birth, T mutation) : _birth(birth), _mutation(mutation) {}
            void operator()(container_linear<L, T>& s) {
                _tmp = std::move(s);
                for(std::size_t i = 0; i != L - 1; ++i) {
                    s[i] = static_cast<T>(1) - _birth[i] + _birth[i] * intrep::pow((static_cast<T>(1) - _mutation) * _tmp[i] - _mutation * _tmp[i+1], 2);
                }
                s.back() = static_cast<T>(1) - _birth.back() + _birth.back() * intrep::pow(_tmp.back(), 2);
            }
    };

    template<std::size_t L, class T>
    void pgf_binary_splitting(container_linear<L, T>& s, container_linear<L, T>& tmp, const container_linear<L, T>& birth, const T& mutation) {
        tmp = std::move(s);
        for(std::size_t i = 0; i != L - 1; ++i) {
            s[i] = static_cast<T>(1) - birth[i] + birth[i] * intrep::pow((static_cast<T>(1) - mutation) * tmp[i] - mutation * tmp[i+1], 2);
        }
        s.back() = static_cast<T>(1) - birth.back() + birth.back() * intrep::pow(tmp.back(), 2);
    }

    template<std::size_t L, class T>
    void pgf_poisson(container_linear<L, T>& s, container_linear<L, T>& tmp, const container_linear<L, T>& fitness, const T& mutation) {
        tmp = std::move(s);
        for(std::size_t i = 0; i != L -1; ++i) {
            s[i] = (1 - mutation) * tmp[i] + mutation * tmp[i+1];
            s[i] = intrep::exp(fitness[i] * (s[i] - 1));
        }
        s.back() = intrep::exp(fitness.back() * (tmp.back() - 1));
    }

    template<std::size_t L, class T>
    T biggest_difference(container_linear<L, T>& s, container_linear<L, T>& tmp) {
        T diff(0);
        for(std::size_t i = 0; i != L; ++i) {
            T t = intrep::abs(s[i] - tmp[i]);
            if(diff < t) {
                diff = t;
            }
        }
        return diff;
    }

    template<std::size_t L, class T>
    T extinction_probability(const container_linear<L, T>& birth, const T& mutation) {
        container_linear<L, T> s;
        container_linear<L, T> tmp;
        s.fill(0);
        std::size_t count = 0;
        do{
            pgf_binary_splitting(s, tmp, birth, mutation);
            ++count;
        } while(count < L * L || biggest_difference(s, tmp) > 1e-300);
        return s.front();
    }

    template<class T,
        class Function,
        class FitnessValues>
    T generic_extinction_probability(Function&& func, FitnessValues&& fit, T&& mutation) {
        constexpr std::size_t L = fit.size();
        container_linear<L, T> s;
        container_linear<L, T> tmp;
        s.fill(0);
        std::size_t count = 0;
        do{
            func(std::ref(s), std::ref(tmp), std::forward<FitnessValues>(fit), std::forward<T>(mutation));
            ++count;
        } while(count < L * L || biggest_difference(s, tmp) > 1e-300);
        return s.front();
    }
}
