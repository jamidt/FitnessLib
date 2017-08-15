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

#include "../container.hpp"
#include "../landscape/generic_landscape.hpp"
#include <array>
#include <random>

namespace intrep {
/* Following "A review on MC simulation methods [...]" Mode et. al.
 * 2008, Math. Biosci 211
 */
/*
 * This class is not supposed to be used anymore!!
template <std::size_t N, class IntType = long>
class multinomial_dist {
  private:
    std::array<IntType, N> _sample;
    IntType _trials;
    const std::array<double, N>& _prob;

  public:
    // Constructor with pointer to probabilitiy array
    multinomial_dist(IntType trials, const std::array<double, N>& prob)
        : _trials(trials), _prob(prob) {}
    // Return array with random numbers from multionomial sampling
    template <class Generator>
    std::array<IntType, N> operator()(Generator& gen) {
        IntType trials_next = _trials;
        double prob_used = 1;
        for (std::size_t i = 0; i != N; ++i) {
            std::binomial_distribution<IntType> binom_dist(
                trials_next, _prob[i] / prob_used);
            _sample[i] = binom_dist(gen);
            prob_used -= _prob[i];
            trials_next -= _sample[i];
        }
        return _sample;
    }
};
*/

/**@ingroup random
 * @brief Multinomial sampling from a absolute distribution
 *
 * Generate the next population via multinomial sampling
 * The Input parameters are the output for the sample of the new
 * generation and the frequencies are optained via using the current
 * generation. Note that this class calculates the next sample from the
 * absolute values of
 */
template <std::size_t L, class IntType,
          template <std::size_t, class> class Container = CubeValues>
class MultinomialSamplingFromAbsolute {
  private:
    std::mt19937_64& _rng_generator;
    const IntType _trials;

  public:
    /**@brief Constructor
     * Constructor with ref to RNG and the number of trials (i.e. the size of
     * the next generation after sampling).
     * @param[in] rng_generator reference to the rng
     * @param[in] trials number of trials.
     */
    MultinomialSamplingFromAbsolute(std::mt19937_64& rng_generator,
                                    IntType trials)
        : _rng_generator(rng_generator), _trials(trials) {}
    /**@brief Get next sample
     * Generate the next generation for the Wright-Fisher dynamics. Note
     * @param[in,out] inout Reference to the population container.
     * This population container is used for generating the new sample and used
     * for storing the new sample.
     */
    // Sampling the next generation. population contains the last
    // generation and out the next one.
    void getSample(Container<L, IntType>& inout) {
        IntType trials_next = _trials;
        double in_size = static_cast<double>(
            std::accumulate(inout.cbegin(), inout.cend(), 0.0));
        for (std::size_t i = 0; i != inout.size(); ++i) {
            std::binomial_distribution<IntType> binom_dist(trials_next,
                                                           inout[i] / in_size);
            in_size -= inout[i];
            inout[i] = binom_dist(_rng_generator);
            trials_next -= inout[i];
        }
    }
};

/**@ingroup random
 * @brief Multinomial sampling from a given fitness landscape
 *
 * The sample has the total population size _population_size, while the
 * input is not bound to any restriction.
 */
template <std::size_t L, class IntType, class RealType,
          template <std::size_t, class> class Container = CubeValues>
class MultinomialSamplingFromFitness {
  private:
    const GenericLandscape<L, RealType, Container>& _fitness;
    std::mt19937_64& _rng_generator;
    const IntType _population_size;

  public:
    /**@brief Constructor
     *
     * This constructs the multinomial sampling class. Note that the
     * number of trials is given in terms of the population size after
     * the sampling process.
     * @param[in] fitness Reference to the fitness container
     * @param[in] rng_generator Reference to the random number generator
     * @param[in] population_size Total population size
     */
    MultinomialSamplingFromFitness(
        const GenericLandscape<L, RealType, Container>& fitness,
        std::mt19937_64& rng_generator, IntType population_size)
        : _fitness(fitness),
          _rng_generator(rng_generator),
          _population_size(population_size) {}
    MultinomialSamplingFromFitness() = delete;

    /**@brief Samping next generation
     *
     * Calculating the population in the next generation via multinomial
     * sampling.
     * @param[in,out] inout Ref to container containing the previous
     * population. The population after samping is stored in this
     * container. It has the total population size, previously set.
     */
    void getSample(Container<L, IntType>& inout) {
        RealType total_fitness(0);
        for (std::size_t vertex = 0; vertex != _fitness.size(); ++vertex) {
            total_fitness += _fitness[vertex] * inout[vertex];
        }
        IntType remaining_population = _population_size;
        RealType prob_unused = static_cast<RealType>(1);
        for (std::size_t vertex = 0; vertex != inout.size(); ++vertex) {
            RealType prob = _fitness[vertex] * inout[vertex] / total_fitness;
            // Make sure, that the probability used is non negative
            if (remaining_population != 0 && prob_unused > static_cast<RealType>(0)) {
                std::binomial_distribution<IntType> binom_dist(
                    remaining_population, prob / prob_unused);
                inout[vertex] = binom_dist(_rng_generator);
                prob_unused -= prob;
                remaining_population -= inout[vertex];
            } else {
                inout[vertex] = 0;
            }
        }
    }
};

}  // End namespace intrep
