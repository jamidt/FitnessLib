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

#include "boost/functional/hash.hpp"
#include "wright_fisher.hpp"
#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <iostream>

namespace intrep {

// Line fo descent path
using LineOfDescent_t = std::vector<std::size_t>;

// Save how often a given pathway is used
template <class Type>
using PathCount = std::tuple<LineOfDescent_t, Type>;

// Vector containing the pathways and the number of times they are taken
template <class Type>
using PathVec = std::vector<PathCount<Type>>;

template <class IntType>
using PathCount_t =
    std::unordered_map<LineOfDescent_t, IntType, boost::hash<LineOfDescent_t>>;

// Container for storing the hashed pathways and their values
template <class T>
using PathHashTable_t = std::vector<std::pair<std::size_t, T>>;

// Search for the value t in vector<T>. Returns true if found
template <class T>
bool inVec(const std::vector<T>& vec, const T& t) {
    for (const auto& v : vec) {
        if (t == v) {
            return true;
        }
    }
    return false;
}

/**
 * TODO: LinesOfDescentHash and LinesOfDescentPath are more or less the
 * same. The only thing that differes is the step, where the pathways
 * are stored. Apart from this, all the other steps should be done in
 * one class!
 */

/**
 * Store the lines of descent with hash value of the pathways.
 */
template <std::size_t L, class IntType = long long, class RealType = double,
          template <std::size_t, class> class Container = CubeValues>
class LinesOfDescentHash {
  private:
    WrightFrisherDynamicsTrack<L, IntType, RealType, Container> _wf_dyn;
    IntType _runs;
    using PathTable_t = PathHashTable_t<IntType>;
    PathTable_t _pathways;
    IntType _max_gen;
    boost::hash<LineOfDescent_t> hash_lod;
    const GenericLandscape<L, RealType, Container>& _fitness;

    // Find line of descent for a given fitness landscape by
    // starting at vertex going back to the wild type
    void findPath(std::size_t vertex) {
        LineOfDescent_t path;
        path.emplace_back(vertex);
        IntType generation = _max_gen;
        const auto& spawn_events = _wf_dyn.getSpawnEvents();
        const auto& mutation_events = _wf_dyn.getMutationEvents();
        // TODO: Put this part in a member function of the
        // WrightFisherDynamicsTrack class that returns the appropriate
        // pathway as vector
        while (vertex != 0) {
            // Note that by definition spawn[vertex] cannot be
            // empty.
            // Search for the last time before the current
            // generation that the site spawned into existence.
            for (const auto& gen : spawn_events.at(vertex)) {
                if (gen < generation) {
                    generation = gen;
                    break;
                } else {
                    continue;
                }
            }
            // Search for the neighbour_vertex from which the
            // mutation happened when the vertex spawned into
            // existence. If this happened from two different
            // vertices at the same generation, this is only attributed
            // to one of them. However, this should be a neglectable
            // incident.
            for (std::size_t k = 0; k != L; ++k) {
                std::size_t neighbour_vertex = vertex ^ (1 << k);
                bool found = inVec(mutation_events.at(neighbour_vertex).at(k),
                                   generation);
                if (found == true) {
                    path.emplace_back(neighbour_vertex);
                    vertex = neighbour_vertex;
                    break;
                }
            }
        }
        // Check if path is already in the pathway vector, if yes
        // add one to the corresponding path, otherwise append path.
        std::size_t hash_path = hash_lod(path);
        auto iter_pw = std::find_if(
            _pathways.begin(), _pathways.end(),
            [hash_path](const auto& x) { return x.first == hash_path; });
        if (iter_pw == _pathways.end()) {
            _pathways.emplace_back(hash_path, 1);
        } else {
            iter_pw->second += 1;
        }
    }

  public:
    template <class S>
    LinesOfDescentHash(const GenericLandscape<L, RealType, Container>& fitness,
                       RealType mutationrate, IntType pop_size, IntType runs,
                       IntType max_gen, S seed)
        : _wf_dyn(fitness, mutationrate, pop_size, seed),
          _runs(runs),
          _max_gen(max_gen) {}
    LinesOfDescentHash(const GenericLandscape<L, RealType, Container>& fitness,
                       RealType mutationrate, IntType pop_size, IntType runs,
                       IntType max_gen)
        : LinesOfDescentHash(fitness, mutationrate, pop_size, runs, max_gen,
                             std::random_device()()) {}

    void findPathways() {
        // Find global maximum of the fitness landscape
        RealType optimum = _fitness.front();
        std::size_t vertex_go = 0;
        for(std::size_t vertex = 0; vertex != _fitness.size(); ++vertex) {
            if(optimum < _fitness[vertex]) {
                optimum = _fitness[vertex];
                vertex_go = vertex;
                continue;
            }
            else {
                continue;
            }
        }
        for (IntType run = 0; run != _runs; ++run) {
            _wf_dyn.reset();
            // Perform the time steps until max_gen is reached or the
            // population fixates (note the population is not fixated,
            // but rather a largest population at one point
            for (IntType gen = 0; gen != _max_gen; ++gen) {
                _wf_dyn.nextGeneration();
                // If the population is at the global maximum of the
                // fitness landscape and is fixed there, terminate the
                // dynamics
                const auto& population = _wf_dyn.getPopulation();
                if(population == _wf_dyn.getPopSize()) {
                    break;
                }
            }
            // Check, which vertex has the largest (sub) population
            IntType largest_pop_type = 0;
            std::size_t vertex_largest_subpop = 0;
            for (std::size_t vertex = 0;
                 vertex != _wf_dyn.getPopulation().size(); ++vertex) {
                auto tmp_sub_pop = _wf_dyn.getPopulation()[vertex];
                if (largest_pop_type < tmp_sub_pop) {
                    largest_pop_type = tmp_sub_pop;
                    vertex_largest_subpop = vertex;
                }
            }
            // Find the pathway from the wild type to the vertex with
            // the largest sub population
            findPath(vertex_largest_subpop);
        }
    }

    std::size_t size() const { return _pathways.size(); }
    const PathTable_t& getPathVec() const { return _pathways; }
    // Get access to the underlying WF dynamics
    const WrightFrisherDynamicsTrack<L, IntType, RealType, Container>&
    getWFDynamics() const {
        return _wf_dyn;
    }
};

/**
 * Store the lines of descent with pathway
 */
template <std::size_t L, class IntType = long long, class RealType = double,
          template <std::size_t, class> class Container = CubeValues>
class LinesOfDescentPath {
  private:
    WrightFrisherDynamicsTrack<L, IntType, RealType, Container> _wf_dyn;
    IntType _runs;
    using PathTable_t = PathCount_t<IntType>;
    PathTable_t _pathways;
    IntType _max_gen;

    // Find line of descent for a given fitness landscape by
    // starting at vertex going back to the wild type
    void findPath(std::size_t vertex) {
        LineOfDescent_t path;
        path.emplace_back(vertex);
        IntType generation = _max_gen;
        const auto& spawn_events = _wf_dyn.getSpawnEvents();
        const auto& mutation_events = _wf_dyn.getMutationEvents();
        // TODO: Put this part in a member function of the
        // WrightFisherDynamicsTrack class that returns the appropriate
        // pathway as vector
        while (vertex != 0) {
            // Note that by definition spawn[vertex] cannot be
            // empty.
            // Search for the last time before the current
            // generation that the site spawned into existence.
            for (const auto& gen : spawn_events.at(vertex)) {
                if (gen < generation) {
                    generation = gen;
                    break;
                } else {
                    continue;
                }
            }
            // Search for the neighbour_vertex from which the
            // mutation happened when the vertex spawned into
            // existence. If this happened from two different
            // vertices at the same generation, this is only attributed
            // to one of them. However, this should be a neglectable
            // incident.
            for (std::size_t k = 0; k != L; ++k) {
                std::size_t neighbour_vertex = vertex ^ (1 << k);
                bool found = inVec(mutation_events.at(neighbour_vertex).at(k),
                                   generation);
                if (found == true) {
                    path.emplace_back(neighbour_vertex);
                    vertex = neighbour_vertex;
                    break;
                }
            }
        }
        // Check if path is already in the pathway vector, if yes
        // add one to the corresponding path, otherwise append path.
        auto iter_pw = _pathways.find(path);
        if (iter_pw == _pathways.end()) {
            _pathways[path] = 1;
        } else {
            iter_pw->second += 1;
        }
    }

  public:
    template <class S>
    LinesOfDescentPath(const GenericLandscape<L, RealType, Container>& fitness,
                       RealType mutationrate, IntType pop_size, IntType runs,
                       IntType max_gen, S seed)
        : _wf_dyn(fitness, mutationrate, pop_size, seed),
          _runs(runs),
          _max_gen(max_gen) {}
    LinesOfDescentPath(const GenericLandscape<L, RealType, Container>& fitness,
                       RealType mutationrate, IntType pop_size, IntType runs,
                       IntType max_gen)
        : LinesOfDescentPath(fitness, mutationrate, pop_size, runs, max_gen,
                             std::random_device()()) {}

    void findPathways() {
        for (IntType run = 0; run != _runs; ++run) {
            _wf_dyn.reset();
            // Perform the time steps until max_gen is reached or the
            // population fixates (note the population is not fixated,
            // but rather a largest population at one point
            for (IntType gen = 0; gen != _max_gen; ++gen) {
                _wf_dyn.nextGeneration();
                if(_wf_dyn.getPopulation().back() == (_wf_dyn.getPopSize() * 0.9)) {
                    break;
                }
            }
            // Check, which vertex has the largest (sub) population
            IntType largest_pop_type = 0;
            std::size_t vertex_largest_subpop = 0;
            for (std::size_t vertex = 0;
                 vertex != _wf_dyn.getPopulation().size(); ++vertex) {
                auto tmp_sub_pop = _wf_dyn.getPopulation()[vertex];
                if (largest_pop_type < tmp_sub_pop) {
                    largest_pop_type = tmp_sub_pop;
                    vertex_largest_subpop = vertex;
                }
            }
            // Find the pathway from the wild type to the vertex with
            // the largest sub population
            findPath(vertex_largest_subpop);
        }
    }

    std::size_t size() const { return _pathways.size(); }
    const PathTable_t& getPathVec() const { return _pathways; }
    // Get access to the underlying WF dynamics
    const WrightFrisherDynamicsTrack<L, IntType, RealType, Container>&
    getWFDynamics() const {
        return _wf_dyn;
    }
};

/**
 * Store the lines of descent with pathway
 */
template <std::size_t L, class IntType = long long, class RealType = double,
          template <std::size_t, class> class Container = CubeValues>
class LinesOfDescentGO {
  private:
    WrightFrisherDynamicsTrack<L, IntType, RealType, Container> _wf_dyn;
    IntType _runs;
    using PathTable_t = PathCount_t<IntType>;
    PathTable_t _pathways;
    using PathDistribution_t = std::vector<RealType>;
    PathDistribution_t _path_dist;
    IntType _max_gen;
    IntType _defect_runs; // Count the number of runs, where the simulation didn't end at the L mutant

    // Find line of descent for a given fitness landscape by
    // starting at vertex going back to the wild type
    void findPath(std::size_t vertex) {
        LineOfDescent_t path;
        path.emplace_back(vertex);
        IntType generation = _max_gen;
        const auto& spawn_events = _wf_dyn.getSpawnEvents();
        const auto& mutation_events = _wf_dyn.getMutationEvents();
        // TODO: Put this part in a member function of the
        // WrightFisherDynamicsTrack class that returns the appropriate
        // pathway as vector
        while (vertex != 0) {
            // Note that by definition spawn[vertex] cannot be
            // empty.
            // Search for the last time before the current
            // generation that the site spawned into existence.
            for (const auto& gen : spawn_events.at(vertex)) {
                if (gen < generation) {
                    generation = gen;
                    break;
                } else {
                    continue;
                }
            }
            // Search for the neighbour_vertex from which the
            // mutation happened when the vertex spawned into
            // existence. If this happened from two different
            // vertices at the same generation, this is only attributed
            // to one of them. However, this should be a neglectable
            // incident.
            for (std::size_t k = 0; k != L; ++k) {
                std::size_t neighbour_vertex = vertex ^ (1 << k);
                bool found = inVec(mutation_events.at(neighbour_vertex).at(k),
                                   generation);
                if (found == true) {
                    path.emplace_back(neighbour_vertex);
                    vertex = neighbour_vertex;
                    break;
                }
            }
        }
        // Check if path is already in the pathway vector, if yes
        // add one to the corresponding path, otherwise append path.
        auto iter_pw = _pathways.find(path);
        if (iter_pw == _pathways.end()) {
            _pathways[path] = 1;
        } else {
            iter_pw->second += 1;
        }
    }

  public:
    template <class S>
    LinesOfDescentGO(const GenericLandscape<L, RealType, Container>& fitness,
                       RealType mutationrate, IntType pop_size, IntType runs,
                       IntType max_gen, S seed)
        : _wf_dyn(fitness, mutationrate, pop_size, seed),
          _runs(runs),
          _max_gen(max_gen), _defect_runs(0) {}
    LinesOfDescentGO(const GenericLandscape<L, RealType, Container>& fitness,
                       RealType mutationrate, IntType pop_size, IntType runs,
                       IntType max_gen)
        : LinesOfDescentGO(fitness, mutationrate, pop_size, runs, max_gen,
                             std::random_device()()) {}

    void findPathways() {
        _path_dist.clear();
        _pathways.clear();
        for (IntType run = 0; run != _runs; ++run) {
            _wf_dyn.reset();
            // Perform the time steps until max_gen is reached or the
            // population fixates at the global optimum (L mutant)
            for (IntType gen = 0; gen != _max_gen; ++gen) {
                _wf_dyn.nextGeneration();
                if(_wf_dyn.getPopulation().back() == _wf_dyn.getPopSize()) {
                    break;
                }
            }
            // Check, which vertex has the largest (sub) population
            IntType largest_pop_type = 0;
            std::size_t vertex_largest_subpop = 0;
            for (std::size_t vertex = 0;
                 vertex != _wf_dyn.getPopulation().size(); ++vertex) {
                auto tmp_sub_pop = _wf_dyn.getPopulation()[vertex];
                if (largest_pop_type < tmp_sub_pop) {
                    largest_pop_type = tmp_sub_pop;
                    vertex_largest_subpop = vertex;
                }
            }
            std::size_t l_mutant = _wf_dyn.getPopulation().size() - 1;
            if(vertex_largest_subpop != l_mutant) {
                _defect_runs++;
            }
            else {
                // Find the pathway from the wild type to the vertex with
                // the largest sub population
                findPath(l_mutant);
            }
        }

        // Crate probability distribution
        IntType total_success_runs = std::accumulate(
            _pathways.cbegin(), _pathways.cend(), static_cast<IntType>(0),
            [](const auto& x, const auto& val) {
                return x + std::get<1>(val);
            });

        for (const auto& path : _pathways) {
            RealType path_prob = static_cast<RealType>(std::get<1>(path)) /
                                 static_cast<RealType>(total_success_runs);
            _path_dist.emplace_back(path_prob);
        }
    }

    std::size_t size() const { return _pathways.size(); }
    const PathTable_t& getPathVec() const { return _pathways; }
    // Get access to the underlying WF dynamics
    const WrightFrisherDynamicsTrack<L, IntType, RealType, Container>&
    getWFDynamics() const { return _wf_dyn; }

    const PathDistribution_t& getPathDistribution() const { return _path_dist; }
};

// Store all results here. Note that the second vector in the tuple
// stores the probability with which the path in the first vector in the
// tuple is taken.
template <class RealType>
using ResultContainerPath =
    std::vector<std::pair<std::vector<std::size_t>, std::vector<RealType>>>;

template <class RealType>
using ResultContainerPath_t =
    std::unordered_map<LineOfDescent_t, std::vector<RealType>,
                       boost::hash<LineOfDescent_t>>;

// This class is for storing multiple results for different runs
// TODO: This container should only store the hash values of the
// pathways!!!
template <class IntType = long long, class RealType = double>
class Result_t {
  private:
    ResultContainerPath_t<RealType> _results;
    std::size_t _length;

  public:
    Result_t() : _length(0) {}
    void appendResults(const PathCount_t<IntType>& res) {
        IntType total_pathways = std::accumulate(
            res.cbegin(), res.cend(), 0,
            [](const IntType& a, const auto& tuple_path_count) -> IntType {
                return a + std::get<1>(tuple_path_count);
            });
        for (const auto& x : res) {
            appendPair(x, total_pathways);
        }
        // If a pathway was not taken in this realzation append 0
        for (auto& x : _results) {
            auto& results_prob_vec = std::get<1>(x);
            if (results_prob_vec.size() < (_length + 1)) {
                results_prob_vec.emplace_back(0);
            }
        }
        ++_length;
    }
    const ResultContainerPath_t<RealType>& getResults() const {
        return _results;
    }

    using iterator = typename ResultContainerPath_t<RealType>::iterator;
    using const_iterator =
        typename ResultContainerPath_t<RealType>::const_iterator;

    iterator begin() { return _results.begin(); }
    iterator end() { return _results.end(); }
    const_iterator begin() const { return _results.cbegin(); }
    const_iterator end() const { return _results.cend(); }
    const_iterator cbegin() const { return _results.cbegin(); }
    const_iterator cend() const { return _results.cend(); }

  private:
    void appendPair(const typename PathCount_t<IntType>::value_type& count,
                    IntType total_pathways) {
        const auto& count_vector = std::get<0>(count);
        const auto& count_value = std::get<1>(count);
        // Check if path is already stored in container
        // If path is in container _result, append the normalized
        // probabilty using total_pathways.
        // Else emplace path in container with zeros for the previous
        // results.
        auto iter = _results.find(count_vector);
        if (iter != _results.end()) {
            iter->second.emplace_back(static_cast<RealType>(count_value) /
                                      static_cast<RealType>(total_pathways));
        } else {
            std::vector<RealType> vec(_length, 0);
            vec.emplace_back(count_value /
                             static_cast<RealType>(total_pathways));
            _results[count_vector] = vec;
        }
    }
};

}  // end namespace intrep
