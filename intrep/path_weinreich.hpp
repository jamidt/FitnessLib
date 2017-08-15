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

/**
 * Calculating the path probabilities according to the paper by
 * Weinreich et. al. (Darwinian Evolution can follow only very few
 * mutational pathways to fitter proteins.
 */

/**
 * Don't use this anymore, use pathweight_adaptive.hpp instead
 */

#pragma once

#include "landscape.hpp"
#include "pathways.hpp"
#include "utility.hpp"
#include "math.hpp"

#include <unordered_map>
#include "boost/functional/hash.hpp" // <- NOTE: Boost has to be installed. No std::hash for std::array :`(

#include <string>
#include <sstream>
#include <bitset>

namespace intrep {

typedef std::pair<std::size_t, std::size_t> pair_inter;

/**
 * Print pair of loci in the form
 * (int(i), binary(i)) -> (int(j), binary(j))
 */
template<std::size_t loci>
std::string str_pair_bin(pair_inter p) {
    std::stringstream stream;
    stream << intrep::str_vertex_bin<loci>(p.first) << " -> " << intrep::str_vertex_bin<loci>(p.second);
    return stream.str();
}


/**
 * Class saves the selection coefficent, fixation probability and step
 * probability at a vertex.
 */
template<class T>
struct values {
    T       fixation_prob;
    T       sel_coeff;
    T       step_prob;
    values() : fixation_prob(0), sel_coeff(0), step_prob(0) {} // Default constr for unordered_map
    values(T sel, T fix, T step) : fixation_prob(fix), sel_coeff(sel), step_prob(step) {}
};

/**
 * Generic class for calculating the path probabilities according to
 * Weinreich.
 */
template<std::size_t loci, class T = double>
class weinreich_path_prob {
    protected:
        std::unordered_map<pair_inter, values<T>, boost::hash<pair_inter>> _pair_value;
    public:
        // Get values for a pair
        values<T> get(pair_inter p) { return _pair_value.find(p); }

        // Calculate the path probability of the given pathway
        // Note that P must have an iterator
        template<class P>
        T path(const P&);

        /// Iterators
        using iterator          = typename std::unordered_map<pair_inter, values<T>, boost::hash<pair_inter>>::iterator;
        using const_iterator    = typename std::unordered_map<pair_inter, values<T>, boost::hash<pair_inter>>::const_iterator;
        iterator begin() { return _pair_value.begin(); }
        iterator end() {return _pair_value.end(); }
        const_iterator cbegin() { return _pair_value.cbegin(); }
        const_iterator cend() { return _pair_value.cend(); }
};

/**
 * Calculate path probability obeying the 0->1 rule
 * Note that this gives wrong and inconsistent results, since exp(...)>1
 * GO AWAY, DO NOT USE!!!!!!
 */
template<std::size_t loci, class T = double>
class weinreich_path_prob_01_full: public weinreich_path_prob<loci, T> {
    public:
        weinreich_path_prob_01_full(const intrep::generic_landscape<loci, T>&);
};

/**
 * Calculate path probability with each step in the direction of
 * increasing fitness and obeying the 0->1 rule
 */
template<std::size_t loci, class T>
class weinreich_path_prob_01 : public weinreich_path_prob<loci, T> {
    public:
        weinreich_path_prob_01(const intrep::generic_landscape<loci, T>&);
};

/**
 * Calculate the path probability with each step in the direction of
 * increasing fitness and obeying the 0->1 rule with population size
 * dependence
 */
template<std::size_t loci, class T>
class weinreich_path_prob_01_pop : public weinreich_path_prob<loci, T> {
    public:
        template<class N>
        weinreich_path_prob_01_pop(const intrep::generic_landscape<loci, T>&, N);
};

/*
template<std::size_t loci, class T = double>
class neighbourhood01{
    private:
        std::unordered_map<std::size_t, values<T>>      _val;
        std::size_t                                     _vertex;
        const intrep::generic_landscape<loci, T>*       _fitness;
    public:
        neighbourhood01(const intrep::generic_landscape<loci, T>&, std::size_t);
        T fix(std::size_t n) { return _val[n].fixation_prob; }
        T sel(std::size_t n) { return _val[n].sel_coeff; }
        T step_pr(std::size_t n) { return _val[n].step_prob; }
};
*/
/*******
 * Functions
 */

/********
 * fixation_probabilities_01
 */

template<std::size_t loci, class T>
template<class P>
T weinreich_path_prob<loci, T>::path(const P& path) {
    typename P::const_iterator      iter1   = path.cbegin();
    typename P::const_iterator      iter2   = path.cbegin();
    ++iter2;
    T       path_prob = static_cast<T>(1);
    for(; iter2 != path.cend(); ++iter2, ++iter1) {
        path_prob       = path_prob * _pair_value.at(std::make_pair(*iter1, *iter2)).step_prob;
    }
    return path_prob;
}

template<std::size_t loci, class T>
weinreich_path_prob_01_full<loci, T>::weinreich_path_prob_01_full(const intrep::generic_landscape<loci, T>& fitness) {
    for(std::size_t vertex = 0; vertex < intrep::verticies(loci); ++vertex) {
        T                       step_fix    = static_cast<T>(0);
        for(std::size_t l = 0; l < loci; l++) {
            std::size_t             step        = 1 << l;
            std::size_t             neighbour   = vertex | step;
            // Iterate through neighbours and calculate
            // the selection coefficient, fixation probability
            // and step probability
            if(neighbour != vertex) {
                pair_inter  ij              = std::make_pair(vertex, neighbour);
                T           selection_coeff = fitness.at(neighbour) / fitness.at(vertex) - static_cast<T>(1);
                T           pr_fixation     = static_cast<T>(1) - intrep::exp(-static_cast<T>(2) * selection_coeff);
                values<T>   val_tmp(selection_coeff, pr_fixation, static_cast<T>(0));
                this->_pair_value.insert(std::make_pair(ij, val_tmp));
                step_fix                    = step_fix + pr_fixation;
            }
        }
        // Calculate the step probability
        for(std::size_t l = 0; l < loci; l++) {
            std::size_t             step        = 1 << l;
            std::size_t             neighbour   = vertex | step;
            if(neighbour != vertex) {
                pair_inter      ij  = std::make_pair(vertex, neighbour);
                values<T>*      v   = &(this->_pair_value[ij]);
                v->step_prob        = v->fixation_prob / step_fix;
            }
        }
    }
}

template<std::size_t loci, class T>
weinreich_path_prob_01<loci, T>::weinreich_path_prob_01(const intrep::generic_landscape<loci, T>& fitness) {
    for(std::size_t vertex = 0; vertex < intrep::verticies(loci); vertex++) {
        T                       step_fix    = static_cast<T>(0);
        for(std::size_t l = 0; l < loci; l++) {
            std::size_t             step        = 1 << l;
            std::size_t             neighbour   = vertex | step;
            // Iterate through neighbours and calculate
            // the selection coefficient, fixation probability
            // and step probability
            if(neighbour != vertex) {
                pair_inter  ij              = std::make_pair(vertex, neighbour);
                T           selection_coeff = fitness.at(neighbour) / fitness.at(vertex) - static_cast<T>(1);
                T           pr_fixation     = static_cast<T>(0);
                if(selection_coeff > 0) {
                    pr_fixation     = static_cast<T>(1) - intrep::exp(-static_cast<T>(2) * selection_coeff);
                }
                values<T>   val_tmp(selection_coeff, pr_fixation, static_cast<T>(0));
                this->_pair_value.insert(std::make_pair(ij, val_tmp));
                step_fix                    = step_fix + pr_fixation;
            }
        }
        // Calculate the step probability
        for(std::size_t l = 0; l < loci; l++) {
            std::size_t             step        = 1 << l;
            std::size_t             neighbour   = vertex | step;
            if(neighbour != vertex) {
                pair_inter      ij  = std::make_pair(vertex, neighbour);
                values<T>*      v   = &(this->_pair_value[ij]);
                // If selection coeff > 0 calculate the step
                // probability, otherwise it is zero
                if(v->sel_coeff > 0) {
                    v->step_prob        = v->fixation_prob / step_fix;
                }
            }
        }
    }
}

template<std::size_t loci, class T>
template<class N>
weinreich_path_prob_01_pop<loci, T>::weinreich_path_prob_01_pop(const intrep::generic_landscape<loci, T>& fitness, N n) {
    for(std::size_t vertex = 0; vertex < intrep::verticies(loci); vertex++) {
        T                       step_fix    = static_cast<T>(0);
        for(std::size_t l = 0; l < loci; l++) {
            std::size_t             step        = 1 << l;
            std::size_t             neighbour   = vertex | step;
            // Iterate through neighbours and calculate
            // the selection coefficient, fixation probability
            // and step probability
            if(neighbour != vertex) {
                pair_inter  ij              = std::make_pair(vertex, neighbour);
                T           selection_coeff = fitness.at(neighbour) / fitness.at(vertex) - static_cast<T>(1);
                T           pr_fixation     = static_cast<T>(0);
                if(selection_coeff > 0) {
                    pr_fixation     = (1 - intrep::exp(-2 * selection_coeff)) / (1 - intrep::exp(-4 * selection_coeff * static_cast<T>(n))) ;
                }
                values<T>   val_tmp(selection_coeff, pr_fixation, static_cast<T>(0));
                this->_pair_value.insert(std::make_pair(ij, val_tmp));
                step_fix                    = step_fix + pr_fixation;
            }
        }
        // Calculate the step probability
        for(std::size_t l = 0; l < loci; l++) {
            std::size_t             step        = 1 << l;
            std::size_t             neighbour   = vertex | step;
            if(neighbour != vertex) {
                pair_inter      ij  = std::make_pair(vertex, neighbour);
                values<T>*      v   = &(this->_pair_value[ij]);
                // If selection coeff > 0 calculate the step
                // probability, otherwise it is zero
                if(v->sel_coeff > 0) {
                    v->step_prob        = v->fixation_prob / step_fix;
                }
            }
        }
    }
}

/*****
 * neighbourhood01
 */
/*
template<std::size_t loci, class T>
neighbourhood01<loci, T>::neighbourhood01(const intrep::generic_landscape<loci, T>& fitness, std::size_t vert) : _vertex(vert), _fitness(&fitness) {
    for(std::size_t l = 0; l < loci; l++) {
        std::size_t     neighbour   = _vertex | (1 << l);
        if(neighbour != _vertex) {
            T       sc      = fitness->at(neighbour) / _fitness->at(_vertex) - static_cast<T>(1);
            T       pr_fix  = static_cast<T>(2) - intrep::exp(- static_cast<T>(2) * sc);
            _val.emplace(std::make_pair(neighbour, values<T>(sc, pr_fix, static_cast<T>(0))));
        }
    }
    T   cpr     = 0;
    for(const auto& v: _val) {
        cpr     = cpr + v.second.fixation_prob;
    }
    for(const auto& v: _val) {
        v.second.step_prob  = v.second.fixation_prob / cpr;
    }
}
*/
} // End namespace intrep
