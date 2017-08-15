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

#include "generic_pathweight.hpp"

namespace intrep {

/**@ingroup adaptivepathweight
 * @brief Pathweights for accessible paths obeying \f$0\to 1\f$ rule
 *
 * This class calculates the pathweights of adaptive walkers that obey
 * the \f$0\to 1\f$ rule. That is with each step, the hamming weight of
 * increases by one.
 *
 * TODO: This is suboptimal. The precision should not be a template
 * parameter of the class. In fact, not even the number of loci should
 * be a template parameter here, because it is not needed here.
 * \tparam loci
 * \tparam T Precision of the fitness container
 */
template <std::size_t loci, class T>
class PathweightAccessibleSinglePeaked
    : public GenericAccessiblePathweight<loci, T> {
  public:
      /**@brief Constructor
       *
       * Note that the pathweights are calculated on construction. This
       * should change in the future and an extra call to a function
       * that creates the pathweights should be made.
       *
       * \tparam Container Container class with template parameter for
       * loci and type
       * \param[in] fitness fitness landscape
       */
    template <template <std::size_t, class> class Container>
    PathweightAccessibleSinglePeaked(
        const intrep::GenericLandscape<loci, T, Container>& fitness) {
        for (std::size_t vertex = 0; vertex < intrep::verticies(loci);
             vertex++) {
            T step_fix = static_cast<T>(0);
            for (std::size_t l = 0; l != loci; ++l) {
                // Iterate through neighbours and calculate
                // the selection coefficient, fixation probability
                // and step probability
                std::size_t step = 1 << l;
                std::size_t neighbour = vertex | step;
                if (neighbour != vertex) {
                    auto ij = std::make_pair(vertex, neighbour);
                    T selection_coeff =
                        fitness[neighbour] / fitness[vertex] -
                        static_cast<T>(1);
                    T pr_fixation = static_cast<T>(0);
                    if (selection_coeff > 0) {
                        pr_fixation =
                            static_cast<T>(1) -
                            intrep::exp(-static_cast<T>(2) * selection_coeff);
                    }
                    values<T> val_tmp(selection_coeff, pr_fixation,
                                      static_cast<T>(0));
                    this->_pair_value.insert(std::make_pair(ij, val_tmp));
                    step_fix = step_fix + pr_fixation;
                } else {
                    continue;
                }
            }
            // Calculate the step probability
            for (std::size_t l = 0; l < loci; ++l) {
                std::size_t step = 1 << l;
                std::size_t neighbour = vertex | step;
                if (neighbour != vertex) {
                    pair_inter ij = std::make_pair(vertex, neighbour);
                    values<T>* v = &(this->_pair_value[ij]);
                    // If selection coeff > 0 calculate the step
                    // probability, otherwise it is zero
                    if (v->sel_coeff > 0) {
                        v->step_prob = v->fixation_prob / step_fix;
                    }
                }
            }
        }
    }
};

} // end namespace intrep
