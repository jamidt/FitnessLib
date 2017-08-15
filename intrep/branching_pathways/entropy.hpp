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
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

#include <bitset>

namespace intrep {

/**
 * Calculating the branching pathweights from the weight type to the
 * escape L mutant.
 *
 * Calculate the probability distribution, entropy and Kullback Leiber
 * Divergence.
 */
template <std::size_t L, class PrecCalcEscape,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class BranchingPathDistSingle {
  public:
    using Result_t = std::vector<double>;

  private:
    const GenericLandscape<L, double>* _fitness;
    double _mutationrate;
    Result_t _path_probabilities;

  public:
    BranchingPathDistSingle(const GenericLandscape<L, double>* fitness,
                            double mutationrate)
        : _fitness(fitness), _mutationrate(mutationrate) {}
    void createDist() {
        _path_probabilities.clear();
        for (std::size_t k = 0; k != L; ++k) {
            std::size_t vertex = 1 << k;
            recursivePathSearch(vertex, 1);
        }
        // Normalization
        double normalization =
            std::accumulate(_path_probabilities.cbegin(),
                            _path_probabilities.cend(), static_cast<double>(0));
        std::for_each(_path_probabilities.begin(), _path_probabilities.end(),
                      [=](double& x) { x = x / normalization; });
    }

    Result_t& getPathDist() const { return _path_probabilities; }
    using iterator = typename Result_t::iterator;
    using const_iterator = typename Result_t::const_iterator;
    iterator begin() { return _path_probabilities.begin(); }
    const_iterator begin() const { return _path_probabilities.cbegin(); }
    const_iterator cbegin() const { return _path_probabilities.cbegin(); }
    iterator end() { return _path_probabilities.end(); }
    const_iterator end() const { return _path_probabilities.cend(); }
    const_iterator cend() const { return _path_probabilities.cend(); }

    double calcPathEntropy() const {
        double entropy = std::accumulate(
            _path_probabilities.cbegin(), _path_probabilities.cend(),
            static_cast<double>(0),
            [](double& ent, const double& p) { return ent - p * std::log(p); });
        return entropy;
    }
    double calcKLDivergence() const {
        double uniform = 1 / static_cast<double>(_fitness->size());
        double kldiv = std::accumulate(
            _path_probabilities.cbegin(), _path_probabilities.cend(),
            static_cast<double>(0), [=](double& ent, const double& p) {
                return ent + p * std::log(p / uniform);
            });
        return kldiv;
    }

  private:
    double step_weight(std::size_t vertex) {
        double weight =
            (*_fitness)[vertex] /
            (static_cast<double>(1) -
             (1 - _mutationrate * hammingweight(vertex)) * (*_fitness)[vertex]);
        return weight;
    }
    void recursivePathSearch(std::size_t vertex, double weight) {
        if (hammingweight(vertex) == L) {
            _path_probabilities.push_back(weight);
        } else {
            for (std::size_t k = 0; k != L; ++k) {
                std::size_t tmp_vertex = vertex | (1 << k);
                if (tmp_vertex != vertex) {
                    double tmp_weight = weight * step_weight(vertex);
                    recursivePathSearch(tmp_vertex, tmp_weight);
                } else {
                    continue;
                }
            }
        }
    }
};

template <class FitnessContainer>
class BranchingPathDistMulti {
  public:
      // Key = endpoint. Value = path probabilies to this endpoint
    using Result_t = std::unordered_map<std::size_t, std::vector<double>>;
    // Key = endpoint, Value = escape probability for a offspring
    // starting at this vertex
    using EscProb_t = std::unordered_map<std::size_t, double>;

  private:
    Result_t _path_probabilities;
    EscProb_t _escape_probability;
    const FitnessContainer* _fitness;
    double _mutationrate;

  public:
    BranchingPathDistMulti(const FitnessContainer* fitness,
                           double mutationrate)
        : _fitness(fitness), _mutationrate(mutationrate) {}
    void createNonNormalized() {
        _path_probabilities.clear();
        _escape_probability.clear();
        createEscapeProb();
        recursivePathSearch(0, 1.0);
    }
    void createDist() {
        createNonNormalized();
        // Normalization
        double norm_const = static_cast<double>(0);
        for (const auto& iter : _path_probabilities) {
            for (const auto& iter_vec : iter.second) {
                norm_const += iter_vec;
            }
        }
        for (auto& iter : _path_probabilities) {
            for (auto& iter_vec : iter.second) {
                iter_vec = iter_vec / norm_const;
            }
        }
    }

    double calcPathEntropy() const {
        double entropy = 0;
        for (const auto& iter : _path_probabilities) {
            for (const auto& iter_vec : iter.second) {
                entropy -= iter_vec * std::log(iter_vec);
            }
        }
        return entropy;
    }

    double calcEndEntropy() const {
        double entropy = 0;
        for (const auto& iter : _path_probabilities) {
            double prob_endpoint = 0;
            for (const auto& iter_vec : iter.second) {
                prob_endpoint += iter_vec;
            }
            entropy -= prob_endpoint * std::log(prob_endpoint);
        }
        return entropy;
    }

    double calcConditionalEntropy() const {
        double entropy = 0;
        for (const auto& iter : _path_probabilities) {
            double prob_endpoint = 0;
            for (const auto& iter_vec : iter.second) {
                prob_endpoint += iter_vec;
            }
            for (const auto& iter_vec : iter.second) {
                entropy += iter_vec * std::log(prob_endpoint / iter_vec);
            }
        }
        return entropy;
    }

    const EscProb_t& getEscProb() const { return _escape_probability; }
    const Result_t& getPathDist() const { return _path_probabilities; }

    using iterator = typename Result_t::const_iterator;
    using const_iterator = typename Result_t::const_iterator;
    iterator begin() { return _path_probabilities.begin(); }
    const_iterator begin() const { return _path_probabilities.cbegin(); }
    const_iterator cbegin() const { return _path_probabilities.cbegin(); }
    iterator end() { return _path_probabilities.end(); }
    const_iterator end() const { return _path_probabilities.cend(); }
    const_iterator cend() const { return _path_probabilities.cend(); }

  private:
    double step_weight(std::size_t vertex) {
        return (*_fitness)[vertex] /
               (static_cast<double>(1) -
                (1 - _mutationrate * hammingweight(vertex)) *
                    (*_fitness)[vertex]);
    }
    void recursivePathSearch(std::size_t vertex, double weight) {
        if ((*_fitness)[vertex] > 1) {
            double tmp_weight = weight *
                                std::pow(_mutationrate, hammingweight(vertex)) *
                                calcEscapeProb(vertex);
            _path_probabilities[vertex].push_back(tmp_weight);
        } else {
            for (std::size_t k = 0; k != _fitness->loci(); ++k) {
                std::size_t next_vertex = vertex | (1 << k);
                if (next_vertex != vertex) {
                    double tmp_weight = weight * step_weight(vertex);
                    recursivePathSearch(next_vertex, tmp_weight);
                }
            }
        }
    }
    void createEscapeProb() {
        for (std::size_t vertex = 0; vertex != _fitness->size(); ++vertex) {
            if ((*_fitness)[vertex] > 1) {
                _escape_probability[vertex] = calcEscapeProb(vertex);
            } else {
                continue;
            }
        }
    }
    double calcEscapeProb(std::size_t vertex) {
        double q = 0.9;
        for (long long t = 0; t != 10000; ++t) {
            q = std::exp(-(*_fitness)[vertex] * (static_cast<double>(1) - q));
        }
        return static_cast<double>(1) - q;
    }
};

} // end namespace intrep
