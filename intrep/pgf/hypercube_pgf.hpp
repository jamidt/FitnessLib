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
#include "../math.hpp"

#include "../../MP/mpreal.h"

namespace intrep {

/**@ingroup pgf
 * @brief Generic pgf class
 *
 * This class provides an interface for all other probability generating
 * functions on the hypercube. The call operator must be implemented by
 * all deriving classes with the signature
 * \code{cpp}
 *   void operator(CalcContainer<L, S>&)
 * \endcode
 * @tparam L sequence length
 * @tparam S Type for solution container
 * @tparam F Type for fitness container and mutation rate
 * @tparam CalcContainer template for the container that contains all
 * calculated values
 * @tparam FitnessContainer template for the fitness container. Usually
 * = CalcContainer
 */
template <std::size_t L, class S, class F = S,
          template <std::size_t, class> class CalcContainer = CubeValues,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class GenericPgf {
  protected:
    /// Pointer to the fitness landscape
    const intrep::GenericLandscape<L, F, FitnessContainer>* fitness;
    F mutation;

  public:
    using calc_t = S;
    /**@brief Constructor
     * @param[in] _fitness Pointer to fitness landscape
     * @param[in] _mutation Mutationrate
     */
    GenericPgf(const intrep::GenericLandscape<L, F, FitnessContainer>& _fitness,
               F _mutation)
        : fitness(&_fitness), mutation(_mutation) {}
    virtual void operator()(CalcContainer<L, S>&) = 0;
    /// Calculates extionction probability of escape type
    S ext_vert(std::size_t vert);
    virtual void init(CalcContainer<L, S>&) = 0;
};

/**
 * pgf for the case where mutations can only appear with increasing
 * fitness.
 */

/**@ingroup pgf
 * @brief Adaptive pgf class
 *
 * This provides the pgf for the case where the population can only move
 * to the next locus with increasing fitness. These are the steps an
 * population would take that moves along adaptive steps.
 *
 * @tparam L sequence length
 * @tparam S Type for solution container
 * @tparam F Type for fitness container and mutation rate
 * @tparam CalcContainer template for the container that contains all
 * calculated values
 * @tparam FitnessContainer template for the fitness container. Usually
 * = CalcContainer
 */
template <std::size_t L, class S, class F = S,
          template <std::size_t, class> class CalcContainer = CubeValues,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class PgfAdaptive
    : public GenericPgf<L, S, F, CalcContainer, FitnessContainer> {
  private:
    CalcContainer<L, S> s;

  public:
    /**@brief Constructor
     *
     * @param[in] _fitness will be forwarded to the pointer that holds
     * the fitness landscape
     * @param[in] _mutation Mutationrate
     */
    PgfAdaptive(
        const intrep::GenericLandscape<L, F, FitnessContainer>& _fitness,
        F _mutation)
        : GenericPgf<L, S, F, CalcContainer, FitnessContainer>(_fitness,
                                                               _mutation) {}
    /**@brief Calculation
     *
     * @param[in,out] sol Use sol for calculation and store the solution
     * in this container
     */
    void operator()(CalcContainer<L, S>& sol);

    /**@brief Initialize the container
     *
     * Note this is only used for the pathway calculations
     */
    void init(CalcContainer<L, S>& sol);
};

/**@ingroup pgf
 * @brief Steps with \f$0\to 1\f$
 *
 * This pgf is used for the case where only \f$0\to 1\f$ mutations can appear,
 * regardless of the fitness difference in this step.
 *
 * @tparam L sequence length
 * @tparam S Type for solution container
 * @tparam F Type for fitness container and mutation rate
 * @tparam CalcContainer template for the container that contains all
 * calculated values
 * @tparam FitnessContainer template for the fitness container. Usually
 * = CalcContainer
 */
template <std::size_t L, class S, class F = S,
          template <std::size_t, class> class CalcContainer = CubeValues,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class Pgf01 : public GenericPgf<L, S, F, CalcContainer, FitnessContainer> {
  private:
    CalcContainer<L, S> s;
    std::size_t tmp, hamming_weight;

  public:
    /**@brief Constructor
     * @param[in] _fitness Pointer to fitness landscape
     * @param[in] _mutation Mutationrate
     */
    Pgf01(const intrep::GenericLandscape<L, F, FitnessContainer>& _fitness,
          F _mutation)
        : GenericPgf<L, S, F, CalcContainer, FitnessContainer>(_fitness,
                                                               _mutation) {}

    /// Calculate the pgf and store value again in sol
    void operator()(CalcContainer<L, S>& sol);

    /// Exrinction probability of the final mutant. I.e. (1,1,...,1)
    void init(CalcContainer<L, S>& sol);
};

template <std::size_t L, class F,
          template <std::size_t, class> class CalcContainer>
decltype(auto) makePgf01(const GenericLandscape<L, F, CalcContainer>& fitness,
                         F mutationrate) {
    return Pgf01<L, F, F, CalcContainer, CalcContainer>{fitness, mutationrate};
}

/**@ingroup pgf
 * @brief \f$0\to 1\f$ with increasing fitness
 *
 * This pgf is used for the case where offspings can only appear with
 * increasing fitness and \f$0\to 1\f$.
 *
 * @tparam L sequence length
 * @tparam S Type for solution container
 * @tparam F Type for fitness container and mutation rate
 * @tparam CalcContainer template for the container that contains all
 * calculated values
 * @tparam FitnessContainer template for the fitness container. Usually
 * = CalcContainer
 */
template <std::size_t L, class S, class F = S,
          template <std::size_t, class> class CalcContainer = CubeValues,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class PgfAdaptive01
    : public GenericPgf<L, S, F, CalcContainer, FitnessContainer> {
  private:
    CalcContainer<L, S> s;

  public:
    /**@brief Constructor
     *
     * @param[in] f Fitness landscape
     * @param[in] m Mutation probability
     */
    PgfAdaptive01(const intrep::GenericLandscape<L, F, FitnessContainer>& f,
                  F m)
        : GenericPgf<L, S, F, CalcContainer, FitnessContainer>(f, m) {}

    /// Calculate the pgf and store values in sol
    void operator()(CalcContainer<L, S>& sol);

    /// Exrinction probability of the final mutant. I.e. (1,1,...,1)
    void init(CalcContainer<L, S>& sol);
};

/**@brief General case with back mutations
 *
 * This is the most general pgf. Here backmutations are allowed from all
 * types.
 *
 * @tparam L sequence length
 * @tparam S Type for solution container
 * @tparam F Type for fitness container and mutation rate
 * @tparam CalcContainer template for the container that contains all
 * calculated values
 * @tparam FitnessContainer template for the fitness container. Usually
 * = CalcContainer
 */
template <std::size_t L, class S, class F = S,
          template <std::size_t, class> class CalcContainer = CubeValues,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class PgfBackmutations
    : public GenericPgf<L, S, F, CalcContainer, FitnessContainer> {
  private:
    CalcContainer<L, S> s;

  public:
    /**@brief Constructor
     *
     * @param[in] _fitness Fitness landscape
     * @param[in] _mutation Mutation probability
     */
    PgfBackmutations(
        const intrep::GenericLandscape<L, F, FitnessContainer>& _fitness,
        F _mutation)
        : GenericPgf<L, S, F, CalcContainer, FitnessContainer>(_fitness,
                                                               _mutation) {}

    /// Calculate the pgf and store values in sol
    void operator()(CalcContainer<L, S>& sol);

    /// Exrinction probability of the final mutant. I.e. (1,1,...,1)
    void init(CalcContainer<L, S>& sol);
};

/**@brief Backmutations except from supercritical types
 *
 * In this case backmutations are allowed except for the supercritical
 * types, i.e. the types that have a survival probability \f$>0\f$.
 *
 * @tparam L sequence length
 * @tparam S Type for solution container
 * @tparam F Type for fitness container and mutation rate
 * @tparam CalcContainer template for the container that contains all
 * calculated values
 * @tparam FitnessContainer template for the fitness container. Usually
 * = CalcContainer
 */
template <std::size_t L, class S, class F = S,
          template <std::size_t, class> class CalcContainer = CubeValues,
          template <std::size_t, class> class FitnessContainer = CubeValues>
class PgfBackmutationsMax : public GenericPgf<L, S, F> {
  private:
    CalcContainer<L, S> s;

  public:
    /**@brief Constructor
     *
     * @param[in] _fitness Fitness landscape
     * @param[in] _mutation Mutation probability
     */
    PgfBackmutationsMax(
        const intrep::GenericLandscape<L, F, FitnessContainer>& _fitness,
        F _mutation)
        : GenericPgf<L, S, F, CalcContainer, FitnessContainer>(_fitness,
                                                               _mutation) {}

    /// Calculate the pgf and store values in sol
    void operator()(CalcContainer<L, S>& sol);

    /// Exrinction probability of the final mutant. I.e. (1,1,...,1)
    void init(CalcContainer<L, S>& sol);
};

/**********************
 * Functions
 */

// GenericPgf
template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
S GenericPgf<L, S, F, CalcContainer, FitnessContainer>::ext_vert(
    std::size_t vert) {
    // This calculates the extinction probability of the vert
    // if it does not mutate to an other vertex. This has to be
    // checked manually.
    S x1, x2(0);
    do {
        x1 = x2;
        x2 = intrep::func<S>::exp((x1 - static_cast<S>(1)) * (*fitness)[vert]);
    } while (
        x2 - x1 >
        static_cast<S>(
            1e-20));  // TODO: Depending on precision of mpreal. **mp_prec_t**
    return x2;
}

// Pgf01
template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void Pgf01<L, S, F, CalcContainer, FitnessContainer>::init(
    CalcContainer<L, S>& sol) {
    sol.fill(static_cast<S>(1));
    sol.back() = this->ext_vert(sol.size() - 1);
}

template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void Pgf01<L, S, F, CalcContainer, FitnessContainer>::operator()(
    CalcContainer<L, S>& sol) {
    s = std::move(sol);
    for (std::size_t vert = 0; vert < sol.size(); vert++) {
        // Calcualte the hamming weight of vertex
        tmp = vert;
        for (hamming_weight = 0; tmp; tmp = tmp >> 1) {
            hamming_weight += tmp & 1;
        }
        sol[vert] = 0;
        tmp = 1;
        // Add L that obey the 0->1 rule
        do {
            if ((vert | tmp) != vert) {
                sol[vert] += s[vert | tmp];
            }
            tmp = tmp << 1;
        } while (tmp < 1 << L);
        sol[vert] = sol[vert] * this->mutation;
        sol[vert] += (1 - this->mutation * (L - hamming_weight)) * s[vert];
        sol[vert] = func<S>::exp((sol[vert] - 1) * (*this->fitness)[vert]);
    }
}

// PgfAdaptive01
template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfAdaptive01<L, S, F, CalcContainer, FitnessContainer>::operator()(
    CalcContainer<L, S>& sol) {
    s = std::move(sol);
    for (std::size_t vertex = 0; vertex < sol.size(); vertex++) {
        sol[vertex] = static_cast<S>(0);
        // Add L that fulfills the 0->1 and increasing fitness rule
        int count = 0;
        for (std::size_t l = 0; l < L; l++) {
            std::size_t neighbour = vertex | (1 << l);
            if (neighbour != vertex) {
                if ((*this->fitness)[vertex] < (*this->fitness)[neighbour]) {
                    sol[vertex] += s[neighbour];
                    count++;
                }
            }
        }
        sol[vertex] *= this->mutation;
        sol[vertex] +=
            (static_cast<S>(1) - this->mutation * static_cast<S>(count)) *
            s[vertex];
        sol[vertex] = intrep::exp((sol[vertex] - static_cast<S>(1)) *
                                  (*this->fitness)[vertex]);
    }
}

template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfAdaptive01<L, S, F, CalcContainer, FitnessContainer>::init(
    CalcContainer<L, S>& sol) {
    sol.fill(1);
    sol.back() = this->ext_vert(sol.size() - 1);
}

// PgfBackmutations
template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfBackmutations<L, S, F, CalcContainer, FitnessContainer>::init(
    CalcContainer<L, S>& sol) {
    sol.fill(1);
    sol.back() = this->ext_vert(sol.size() - 1);
}

template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfBackmutations<L, S, F, CalcContainer, FitnessContainer>::operator()(
    CalcContainer<L, S>& sol) {
    s = std::move(sol);
    for (std::size_t vert = 0; vert < sol.size(); vert++) {
        sol[vert] = 0;
        for (std::size_t k = 0; k < L; k++) {
            sol[vert] += s[vert ^ (1 << k)];
        }
        sol[vert] *= this->mutation;
        sol[vert] += (1 - L * this->mutation) * s[vert];
        sol[vert] =
            intrep::func<S>::exp((sol[vert] - 1) * (*this->fitness)[vert]);
    }
}

// PgfAdaptive
template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfAdaptive<L, S, F, CalcContainer, FitnessContainer>::init(
    CalcContainer<L, S>& sol) {
    // Function searches for fitness values > 1 and initializes
    // sol with the extinction probability for this type.
    for (std::size_t vert = 0; vert < sol.size(); vert++) {
        if ((*this->fitness)(vert) > 1) {
            sol[vert] = this->ext_vert(vert);
        } else {
            sol[vert] = 1;
        }
    }
}

template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfAdaptive<L, S, F, CalcContainer, FitnessContainer>::operator()(
    CalcContainer<L, S>& sol) {
    s = std::move(sol);
    for (std::size_t vert = 0; vert < sol.size(); vert++) {
        sol[vert] = 0;
        std::size_t neighbours = 0;
        for (std::size_t k = 0; k < L; k++) {
            if ((*this->fitness)(vert) < (*this->fitness)[vert ^ (1 << k)]) {
                sol[vert] += s[vert ^ (1 << k)];
                neighbours++;
            }
        }
        sol[vert] = sol[vert] * this->mutation;
        sol[vert] += (1 - this->mutation * neighbours) * s[vert];
        sol[vert] =
            intrep::func<S>::exp((sol[vert] - 1) * (*this->fitness)[vert]);
    }
}

// PgfBackmutationsMax
template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfBackmutationsMax<L, S, F, CalcContainer, FitnessContainer>::operator()(
    CalcContainer<L, S>& sol) {
    s = std::move(sol);
    for (std::size_t vert = 0; vert < sol.size(); vert++) {
        if ((*this->fitness)[vert] <= static_cast<F>(1)) {
            sol[vert] = static_cast<S>(0);
            for (std::size_t n = 0; n < L; n++) {
                sol[vert] += s[vert ^ (1 << n)];
            }
            sol[vert] *= this->mutation;
            sol[vert] += (1 - L * this->mutation) * s[vert];
            sol[vert] = intrep::exp((sol[vert] - 1) * (*this->fitness)[vert]);
        } else {
            sol[vert] = intrep::exp((*this->fitness)[vert] * (s[vert] - 1));
        }
    }
}

template <std::size_t L, class S, class F,
          template <std::size_t, class> class CalcContainer,
          template <std::size_t, class> class FitnessContainer>
void PgfBackmutationsMax<L, S, F, CalcContainer, FitnessContainer>::init(
    CalcContainer<L, S>& sol) {
    for (std::size_t vert = 0; vert < sol.size(); vert++) {
        if ((*this->fitness)[vert] > static_cast<F>(1)) {
            sol[vert] = this->ext_vert(vert);
        } else {
            sol[vert] = static_cast<S>(1);
        }
    }
}

}  // end namespace intrep
