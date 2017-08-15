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

#include <algorithm>

#include "../math.hpp"
#include "../landscape/generic_landscape.hpp"
#include <iostream>

namespace intrep {

/**
 * Caution: is the escape type in PathContainer? Note that the
 * analytical solution requires that the weight stops before the escape
 * type \tau. I.e. the pathway is (\hat{0}\to\sigma^1\to\dots\to\sigma^k)\tau
 */
template <class FitnessContainer, class PathContainer>
auto branching_pathweight(const FitnessContainer& fitness,
                          const PathContainer& path) ->
    typename FitnessContainer::value_type {
    /*
    std::cout << "Pathweight for path:\n";
    for(const auto& p: path) {
        std::cout << p << " ";
    }
    std::cout << std::endl;
    */
    auto end = path.cend();
    --end;
    typename FitnessContainer::value_type weight = 1;
    for (auto iter = path.cbegin(); iter != end; ++iter) {
        weight = weight * fitness[*iter] / (1 - fitness[*iter]);
    }

    return weight;
}

/**@ingroup
 * @brief Calculate pathweight of given path
 *
 * The PathContainer must contain all verticies, with start point *and*
 * end point.
 *
 * @tparam PathContainer Container with all vertices of the path
 * @tparam FitnessContainer Container with fitness values
 */
template <class FitnessContainer, class PathContainer, class MutationRate>
auto calcBranchingPathweigt(const FitnessContainer& fitness,
                            const PathContainer& path, const MutationRate& m) ->
    typename FitnessContainer::value_type {
    // Don't count the first and last step in the path to to weight
    auto end = path.cend();
    --end;
    auto iter = path.cbegin();
    ++iter;
    typename FitnessContainer::value_type weight = 1;
    for (; iter != end; ++iter) {
        weight = fitness[*iter] /
                 (1 - (1 - m * hammingweight(*iter)) * fitness[*iter]);
    }
    return weight;
}

/**@ingroup
 * @brief Calculate weight ratio of a set of pathways
 *
 */
template <class FitnessContainer, class PathwaysContainer, class MutationRate>
decltype(auto) calcBranchingPathweigtRatio(const FitnessContainer& fitness,
                                 const PathwaysContainer& pathways,
                                 const MutationRate& m) {
//    typename FitnessContainer::value_type {
    using calc_t = typename FitnessContainer::value_type;

    auto path_iter = pathways.cbegin();
    calc_t maximal_weight = calcBranchingPathweigt(fitness, *path_iter, m);
    calc_t minimal_weight = maximal_weight;

    ++path_iter;
    for (; path_iter != pathways.cend(); ++path_iter) {
        calc_t tmp_weight = calcBranchingPathweigt(fitness, *path_iter, m);
        if (tmp_weight < minimal_weight) {
            minimal_weight = tmp_weight;
        } else if (maximal_weight < tmp_weight) {
            maximal_weight = tmp_weight;
        }
    }

    return maximal_weight / minimal_weight;
}

template <std::size_t L, class FitnessValue, class PrecCalc = FitnessValue>
class ApproxPathProbSingleEscape {
  private:
    const GenericLandscape<L, FitnessValue>* const _fitness;
    FitnessValue _m;
    FitnessValue _escape_prob;

  public:
    ApproxPathProbSingleEscape(const GenericLandscape<L, FitnessValue>* fitness,
                               FitnessValue m)
        : _fitness(fitness), _m(m) {}

    template <class PathContainer>
    FitnessValue calcPathProb(const PathContainer* path) {
        auto iter = path->cbegin();
        auto path_iter_end = path->cend();
        --path_iter_end;
        FitnessValue path_prob = static_cast<FitnessValue>(1);
        for (; iter != path_iter_end; ++iter) {
            path_prob *=
                (*_fitness)[*iter] /
                (static_cast<FitnessValue>(1) -
                 (static_cast<FitnessValue>(1) - _m * hammingweight(*iter)) *
                     (*_fitness)[*iter]);
        }
        path_prob *= intrep::func<FitnessValue>::pow(_m, L);
        return path_prob * _escape_prob;
    }

    void calcEscapeTypeProb() {
        PrecCalc s(static_cast<PrecCalc>(0));
        for (int i = 0; i != 1000; ++i) {
            s = intrep::func<PrecCalc>::exp((s - static_cast<PrecCalc>(1)) *
                                            _fitness->back());
        }
        s = static_cast<PrecCalc>(1) - s;
        _escape_prob = static_cast<FitnessValue>(s);
    }
};
}
