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
/**
 * This class calculates the path probability of a given pathway.
 * Note that the number of starting individuals has to be taken into
 * account in the pgf.
 */

#include <fstream>
#include "../container.hpp"
#include "../pathways.hpp"
#include <vector>
#include "../math.hpp"

namespace intrep {

template <std::size_t L, class S,
          template <std::size_t, class> class CalcContainer = CubeValues>
class BranchingTimeDist01 {
  private:
    CalcContainer<L, S> s_t2, s_t1;
    S cumul_path_prob;
    S tmp;
    std::size_t time;
    std::vector<S> tmp_advanced;

  public:
    BranchingTimeDist01()
        : cumul_path_prob(static_cast<S>(0)), time(0) {}

    /**@brief Calculate the escape probability at the next timestep
     *
     * @param[in] pgf Reference to the proability generating function
     * @param[in] pathway Reference to the pathway
     * @return Escape probability at the next timestep
     */
    template <class PGF, class P>
    S nextTimestep(PGF& pgf, const P& pathway) {
        if (time != 0) {
            for (std::size_t i = 0; i < pathway.size(); i++) {
                tmp_advanced.at(i) = s_t1.at(pathway.at(i));
            }
            s_t1 = s_t2;
            for (std::size_t i = 0; i < pathway.size(); i++) {
                s_t1[pathway[i]] = tmp_advanced[i];
            }
            pgf(s_t2);
            pgf(s_t1);
        } else {
            // t = 1 timestep
            initialize(pgf);
            pgf(s_t1);
            tmp_advanced.resize(pathway.size());
        }
        time++;
        tmp = s_t2[0] - s_t1[0];
        cumul_path_prob += tmp;
        return tmp;
    }

    /// Returns the current generation
    std::size_t getTime() { return time; }

    /**@brief Get cumulative probability
     *
     * Calculates the cumulative probability of the given pathway up to
     * the current generation.
     * @return Cumulative probability
     */
    S getCProb() { return cumul_path_prob; }

  private:
    template <class PGF>
    void initialize(PGF& pgf) {
        pgf.init(s_t1);
        s_t2 = s_t1;
    }
};


/**@ingroup branchingpath
 * @brief Calculate the path probability of a given pathway
 *
 * TODO: This is reduntand!!!!
 */
template <std::size_t L, class S,
          template <std::size_t, class> class CalcContainer = CubeValues>
class BranchingSinglePathwayProb01 {
    /*
     * TODO:
     *   - Use iterators for the pathway container instead of random
     *     access.
     */
  private:
    S N;
    CalcContainer<L, S> s_t2, s_t1;
    S cumul_path_prob;
    S tmp;
    std::size_t time;
    std::vector<S> tmp_advanced;

  public:
    /**@brief Constructor
     * @param[in] _N Number of individuals starting the process
     */
    template <class P>
    BranchingSinglePathwayProb01(P _N)
        : N(_N), cumul_path_prob(0), time(0) {}
    /// Default constructor, starting the process with one individual
    BranchingSinglePathwayProb01() : BranchingSinglePathwayProb01(1) {}

    /**@brief Calculate the escape probability at the next timestep
     *
     * @param[in] pgf Reference to the proability generating function
     * @param[in] pathway Reference to the pathway
     * @return Escape probability at the next timestep
     */
    template <class PGF, class P>
    S next_timestep(PGF& pgf, const P& pathway);

    /// Returns the current generation
    std::size_t getTime() { return time; }

    /**@brief Get cumulative probability
     *
     * Calculates the cumulative probability of the given pathway up to
     * the current generation.
     * @return Cumulative probability
     */
    S getCProb() { return cumul_path_prob; }

  private:
    template <class PGF>
    void initialize(PGF&);
};




// Do not use
// template<std::size_t L, class S>
// class generic_multi_path {
// protected:
// intrep::CubeValues<L, S> s_t2, s_t1;
// intrep::GenericPathways<L>* pathw;
// long N;
// std::size_t time;
// std::vector<S> cprob; // Cumulative probability of each pathway
// std::vector<S> tmp_advanced;
// std::vector<S> out;
// public:
//// The initial conditions have to be calculated in the
//// constructor of the derived classes.
// generic_multi_path(intrep::GenericPathways<L>& _pathw, long _N)
//: pathw(&_pathw), N(_N), time(1), cprob(_pathw.size())
//, tmp_advanced(_pathw.longest_path()), out(_pathw.size()) {}

// template<class F>
// std::vector<S> next_timestep(F&);

// std::vector<S> get_CProb() { return cprob; }
//};

/****************
 * Functions
 */
template <std::size_t L, class S,
          template <std::size_t, class> class CalcContainer>
template <class PGF, class P>
S BranchingSinglePathwayProb01<L, S, CalcContainer>::next_timestep(
    PGF& pgf,
    const P& pathway) {  // Calculating the pathway probability at the next timestep
    if (time != 0) {
        for (std::size_t i = 0; i < pathway.size(); i++) {
            tmp_advanced.at(i) = s_t1.at(pathway.at(i));
        }
        s_t1 = s_t2;
        for (std::size_t i = 0; i < pathway.size(); i++) {
            s_t1[pathway[i]] = tmp_advanced[i];
        }
        pgf(s_t2);
        pgf(s_t1);
    } else {
        // t = 1 timestep
        initialize(pgf);
        pgf(s_t1);
        tmp_advanced.resize(pathway.size());
    }
    time++;
    tmp = intrep::func<S>::pow(s_t2[0], N) - intrep::func<S>::pow(s_t1[0], N);
    cumul_path_prob += tmp;
    return tmp;
}

template <std::size_t L, class S,
          template <std::size_t, class> class CalcContainer>
template <class PGF>
void BranchingSinglePathwayProb01<L, S, CalcContainer>::initialize(PGF& pgf) {
    pgf.init(s_t1);
    s_t2 = s_t1;
}

// template<std::size_t L, class S>
// template<class F>
// std::vector<S> generic_multi_path<L, S>::next_timestep(F& pgf) {
// for(std::size_t p = 0; p < pathw->size(); p++) {
// for(std::size_t i = 0; i < pathw->at(p).size(); i++) {
// tmp_advanced.at(i) = s_t1.at(pathw->at(p).at(i));
//}
// s_t1 = s_t2;
// for(std::size_t i = 0; i < pathw->at(p).size(); i++) {
// s_t1.at(pathw->at(p).at(i)) = tmp_advanced.at(i);
//}
// pgf(s_t1);
// out.at(p) = intrep::func<S>::pow(s_t1[0], N);
//}
// pgf(s_t2);
// for(std::size_t p = 0; p < pathw->size(); p++) {
// out.at(p) = intrep::func<S>::pow(s_t2, N) - out.at(p);
// cprob.at(p) += out.at(p);
//}
// time++;
// return std::move(out);
//}
}
