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

#include "../landscape/generic_landscape.hpp"
#include "../math.hpp"
#include "../random.hpp"
#include <array>
#include <random>
#include <vector>

namespace intrep {
/**
 * Generic class for the Wright Fisher process. This class acts only as
 * an interface for specializing.
 */
template <std::size_t L, class IntType = long long, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class GenericWrightFrisherDynamics {
  protected:
    const GenericLandscape<L, T, Container>& _fitness;
    const T _mutationrate;
    const IntType _population_size;
    IntType _generation;
    Container<L, IntType> _population;

  public:
    GenericWrightFrisherDynamics() = delete;
    GenericWrightFrisherDynamics(
        const GenericLandscape<L, T, Container>& fitness, T mutationrate,
        IntType pop_size)
        : _fitness(fitness),
          _mutationrate(mutationrate),
          _population_size(pop_size) {}
    void reset() { _generation = 0; }
    T getTotalFitness() const {
        T total_fitness(0);
        for (std::size_t vertex = 0; vertex != _population.size(); ++vertex) {
            total_fitness += _fitness[vertex] * _population[vertex];
        }
        return total_fitness;
    }
    T getMeanFitness() {
        T mean = getTotalFitness();
        T real_pop_size =
            std::accumulate(_population.cbegin(), _population.cend(), 0.0);
        mean /= real_pop_size;
        return mean;
    }
    IntType getPopSize() const { return _population_size; }
    const Container<L, IntType>& getPopulation() const { return _population; }
    IntType getGeneration() const { return _generation; }

    /*
    bool isFixed(IntType error_threshold = 0) {
        if (_population.back() >= _population_size - error_threshold) {
            return true;
        } else {
            return false;
        }
    }
    */
    IntType getGeneration() { return _generation; }
};

// Keeping track of mutation events.
// a[i] gives an array of mutation events from vertex i to the L
// neighbours. a[vertex] = std::array<std::vector<IntType>, L>
// I.e. a[vertex][j] gives the times when a mutation at vertex at loci j happens
// and
// n(vertex) = 0 at that time. Thus, this provides the times where a population
// at j might spawn.
// Note that the neighbouring vertex is given by calculating vertex ^ (1<<j)
template <std::size_t L, class IntType>
using MutationEvents =
    std::array<std::array<std::vector<IntType>, L>, verticies(L)>;

// Keeping track of spawn events on the sites.
template <std::size_t L, class IntType>
using SpawnEvents = std::array<std::vector<IntType>, verticies(L)>;

/**
 * This class is a Wright Fisher process which can be used to track the
 * spawning and mutation to different types. It is used for calculating
 * the pathway probabilities within the WF process.
 */
template <std::size_t L, class IntType = long long, class T = double,
          template <std::size_t, class> class Container = CubeValues>
class WrightFrisherDynamicsTrack
    : public GenericWrightFrisherDynamics<L, IntType, T, Container> {
  private:
    std::mt19937_64 _rng_gen;
    // Number of mutations present in this generation. Mean = mutationrate * L *
    // popsize
    std::poisson_distribution<IntType> _mut_dist;
    // Uniform int in [0, L), used for the mutation step
    std::uniform_int_distribution<std::size_t> _uniform_int;
    MutationEvents<L, IntType> _mutation_events;
    SpawnEvents<L, IntType> _spawn_events;
    Container<L, IntType> _population_previous_generation;
    MultinomialSamplingFromFitness<L, IntType, T, Container>
        _multinomial_distribution;

  public:
    template <class S>
    WrightFrisherDynamicsTrack(const GenericLandscape<L, T, Container>& fitness,
                               T mutationrate, IntType pop_size, S seed)
        : GenericWrightFrisherDynamics<L, IntType, T, Container>::
              GenericWrightFrisherDynamics(fitness, mutationrate, pop_size),
          _rng_gen(seed),
          _mut_dist(mutationrate * pop_size * L),
          _uniform_int(0, L - 1),
          _multinomial_distribution(fitness, _rng_gen, pop_size) {
        reset();
    }
    WrightFrisherDynamicsTrack(const GenericLandscape<L, T, Container>& fitness,
                               T mutationrate, IntType pop_size)
        : WrightFrisherDynamicsTrack(fitness, mutationrate, pop_size,
                                     std::random_device{}()) {}
    void reset() {
        GenericWrightFrisherDynamics<L, IntType, T, Container>::reset();
        this->_population.fill(0);
        this->_population.front() = this->_population_size;
        for (auto& mtype : _mutation_events) {
            for (auto& mtimes : mtype) {
                mtimes.clear();
            }
        }
        for (auto& stype : _spawn_events) {
            stype.clear();
        }
    }
    void nextGeneration() {
        ++(this->_generation);
        mutation();
        selection();
    }
    const MutationEvents<L, IntType>& getMutationEvents() const {
        return _mutation_events;
    }
    const SpawnEvents<L, IntType>& getSpawnEvents() const {
        return _spawn_events;
    }

  private:
    // Performing the mutation step and storing the "new" population in
    // _population. Note this leaves _population in a non constant
    // population size state. That is, the new population size is the original
    // one plus the total number of mutations obtained in this step. The
    // function void selection() has the task of normalizing the population again via multinomial sampling.
    void mutation() {
        // Keeping track of the population before selection. This is
        // used in this function to keep track of newly spawned
        // (occupied) sites in the population.
        _population_previous_generation = this->_population;
        // Generate the number of mutations that appear in this
        // generation.
        IntType n_mutations = _mut_dist(_rng_gen);
        for (std::size_t vertex = 0; vertex != this->_population.size();
             ++vertex) {
            double vertex_mutations =
                this->_population[vertex] * n_mutations /
                static_cast<double>(this->_population_size);
            // Distribute mutations uniformly to all L neighbour types
            // with hamming distance 1.
            for (double x = 1; x <= vertex_mutations; ++x) {
                std::size_t mutation_loci = _uniform_int(_rng_gen);
                std::size_t neighbour_vertex = vertex ^ (1 << mutation_loci);
                // Note that mutation_loci is the loci at which the
                // vertex and its neighbour differ.
                // If neighbour_vertex was unoccupied before mutation,
                // save generation at which the mutation happens.
                if (this->_population[neighbour_vertex] == 0) {
                    _mutation_events[vertex][mutation_loci].emplace_back(
                        this->_generation);
                }
                this->_population[neighbour_vertex] += 1;
            }
        }
    }
    void selection() {
        // Changing the frequencies of the sites according to their
        // fitness.
        _multinomial_distribution.getSample(this->_population);
        // Keeping track of the newly spawned sites
        for (std::size_t vertex = 0; vertex != this->_population.size();
             ++vertex) {
            if (this->_population[vertex] != 0 &&
                _population_previous_generation[vertex] == 0) {
                _spawn_events[vertex].emplace_back(this->_generation);
            }
        }
        // If mutation appeared to this vertex but the selection
        // process selected against the new mutation, delete the
        // mutational process from the mutation_event list
        for (std::size_t vertex = 0; vertex != this->_population.size();
             ++vertex) {
            for (std::size_t k = 0; k != L; ++k) {
                std::size_t to_vertex = vertex ^ (1 << k);
                if (_mutation_events[vertex][k].size() != 0 &&
                    _mutation_events[vertex][k].back() == this->_generation &&
                    !(_spawn_events[to_vertex].size() != 0 &&
                      _spawn_events[to_vertex].back() == this->_generation)) {
                    _mutation_events[vertex][k].pop_back();
                }
            }
        }
        _population_previous_generation = this->_population;
    }
};

// Create WrightFisherDynamicsTrack class (circumvading the lack of
// template deduction for classes pre C++17)
template <std::size_t L, class IntType, class T,
          template <std::size_t, class> class Container = CubeValues>
WrightFrisherDynamicsTrack<L, IntType, T, Container> make_WF_dynamic_track(
    const Container<L, T>& fitness, T mutationrate, IntType pop_size) {
    return WrightFrisherDynamicsTrack<L, IntType, T, Container>(
        fitness, mutationrate, pop_size);
}
}
