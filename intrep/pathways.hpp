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

#include <vector>
#include "utility.hpp"
#include <fstream>
#include <bitset>
#include "landscape/generic_landscape.hpp"

#include <sstream>
#include <string>

namespace intrep {

using Path_t = std::vector<std::size_t>;

// Return a string of a pathway in binary from
template <std::size_t loci>
std::string strPathBin(const Path_t&);

// Return a string of a pathway in the form (int, binary)
template <std::size_t loci>
std::string strPathIntBin(const Path_t&);

// Return a string of a pathway in integer form
template <size_t loci>
std::string strPathInt(const Path_t&);

/**
 * Generic pathway class for storing all possible pathways.
 * The pathways should be constructed using the constructor or by a
 * function called by the constructor (void createPathways(...)).
 */
template <std::size_t loci>
class GenericPathways {
  protected:
    std::vector<Path_t> pathways;

  public:
    typedef typename std::vector<Path_t>::iterator iterator;
    typedef typename std::vector<Path_t>::const_iterator const_iterator;
    typedef Path_t value_type;
    iterator begin() { return pathways.begin(); }
    iterator end() { return pathways.end(); }
    const_iterator begin() const { return pathways.cbegin(); }
    const_iterator end() const { return pathways.cend(); }
    const_iterator cbegin() const { return pathways.cbegin(); }
    const_iterator cend() const { return pathways.cend(); }

    Path_t& operator[](std::size_t i) noexcept { return pathways[i]; }
    const Path_t& operator[](std::size_t i) const noexcept {
        return pathways[i];
    }
    Path_t& at(std::size_t i) { return pathways.at(i); }
    const Path_t& at(std::size_t i) const { return pathways.at(i); }

    std::size_t longest_path();
    std::size_t size() { return pathways.size(); }

    bool printToFile(const char*);
    bool printToFile(std::string filename) {
        return printToFile(filename.c_str());
    }
    // Prints the pathway to a file in human readable form.
    bool printToFileBin(const char*);
    bool printToFileBin(std::string filename) {
        return printToFileBin(filename.c_str());
    }
};

/**
 * Create pathways that go from the wild type to the L mutant (i.e. the
 * one with all mutations. Note, that the pathways are independent of
 * the fitness. There are exactly L! pathways created.
 */
template <std::size_t loci>
class Pathways01 : public GenericPathways<loci> {
  public:
    Pathways01();

  private:
    void createPathways(std::vector<std::size_t>, std::size_t);
};

/**
 * This class calculates all pathways along increasing fitness values.
 * This gives pathways that end at a local fitness maimum. Note that the
 * path can perform backsteps.
 */
template <std::size_t loci>
class PathwaysInc : public GenericPathways<loci> {
  public:
    PathwaysInc(){};
    template <class F, template <std::size_t, class> class FitnessContainer>
    PathwaysInc(const intrep::GenericLandscape<loci, F, FitnessContainer>& f) {
        createPathways(f);
    }
    template <class F, template <std::size_t, class> class FitnessContainer>
    void createPathways(
        const intrep::GenericLandscape<loci, F, FitnessContainer>&);

  private:
    template <class F, template <std::size_t, class> class FitnessContainer>
    void single_path_finder(
        const intrep::GenericLandscape<loci, F, FitnessContainer>*, Path_t,
        std::size_t);
};

/**
 * Create all pathways that follow the 0->1 rule (i.e. no backmutations)
 * and are ascending in fitness values.
 */
template <std::size_t loci>
class PathwaysInc01 : public GenericPathways<loci> {
  public:
    PathwaysInc01() {}
    template <class F, template <std::size_t, class> class FitnessContainer>
    PathwaysInc01(
        const intrep::GenericLandscape<loci, F, FitnessContainer>& f) {
        createPathways(f);
    }
    template <class F, template <std::size_t, class> class FitnessContainer>
    void createPathways(
        const intrep::GenericLandscape<loci, F, FitnessContainer>&);

  private:
    template <class T>
    void single_path_finder(T*, Path_t, std::size_t);
};

template <size_t L, class F,
          template <std::size_t, class> class FitnessContainer>
auto make_pathways_inc(
    const intrep::GenericLandscape<L, F, FitnessContainer>& fitness)
    -> PathwaysInc<L> {
    return PathwaysInc<L>(fitness);
}

template <size_t L, class F,
          template <std::size_t, class> class FitnessContainer>
auto make_pathways_inc_01(
    const intrep::GenericLandscape<L, F, FitnessContainer>& fitness)
    -> PathwaysInc01<L> {
    return PathwaysInc01<L>(fitness);
}

/******************
 * Functions
 */

template <std::size_t loci>
std::string strPathBin(const Path_t& path) {
    std::stringstream stream;
    auto iter = path.cbegin();
    stream << str_vertex_bin<loci>(*iter);
    iter++;
    for (; iter != path.cend(); iter++) {
        stream << " -> " << str_vertex_bin<loci>(*iter);
    }
    return stream.str();
}

template <std::size_t loci>
std::string strPathIntBin(const Path_t& path) {
    std::stringstream stream;
    auto iter = path.cbegin();
    stream << str_vertex_int_bin<loci>(*iter);
    iter++;
    for (; iter != path.cend(); iter++) {
        stream << " -> " << str_vertex_int_bin<loci>(*iter);
    }
    return stream.str();
}

template <std::size_t loci>
std::string strPathInt(const Path_t& path) {
    std::stringstream stream;
    auto iter = path.cbegin();
    stream << str_vertex_bin<loci>(*iter);
    iter++;
    for (; iter != path.cend(); iter++) {
        stream << " -> " << *iter;
    }
    return stream.str();
}

template <std::size_t loci>
std::size_t GenericPathways<loci>::longest_path() {
    std::size_t length = 0;
    for (const auto& p : pathways) {
        if (length < p.size()) {
            length = p.size();
        }
    }
    return length;
}

template <std::size_t loci>
bool GenericPathways<loci>::printToFile(const char* filename) {
    std::ofstream pathways_data(filename);
    pathways_data << "# Pathway number: pathway..." << std::endl;
    for (std::size_t n = 0; n < pathways.size(); n++) {
        pathways_data << n << ": ";
        for (const auto& p : pathways[n]) {
            pathways_data << p << " ";
        }
        pathways_data << "\n";
    }
    pathways_data.close();
    return true;
}

template <std::size_t loci>
bool GenericPathways<loci>::printToFileBin(const char* filename) {
    std::ofstream data_pathways(filename);
    data_pathways << "# Pathway number: pathway..." << std::endl;
    for (std::size_t n = 0; n < pathways.size(); n++) {
        data_pathways << n << ": ";
        for (const auto& v : pathways[n]) {
            data_pathways << std::bitset<loci>(v) << " ";
        }
        data_pathways << std::endl;
    }
    data_pathways.close();
    return true;
}

template <std::size_t loci>
Pathways01<loci>::Pathways01() {
    std::vector<std::size_t> path(loci + 1);
    path[0] = 0;
    createPathways(path, 1);
    GenericPathways<loci>::pathways.shrink_to_fit();
}

template <std::size_t loci>
void Pathways01<loci>::createPathways(std::vector<std::size_t> path,
                                      std::size_t step) {
    // At the end, set (1,1,...,1) = ct::verticies(loci) as the end point of the
    // path
    // and add calculated path to the pathways.
    if (step == loci) {
        path[loci] = intrep::verticies(loci) - 1;
        path.shrink_to_fit();
        GenericPathways<loci>::pathways.emplace_back(path);
    } else {
        for (std::size_t n = 0; n != loci; n++) {
            // Test if the next step obeys the 0->1 mutation rule. If it is
            // fulfilled
            // take step.
            if ((path[step - 1] | 1 << n) != path[step - 1]) {
                path[step] = path[step - 1] | 1 << n;
                createPathways(path, step + 1);
            }
        }
    }
}

template <std::size_t loci>
template <class F, template <std::size_t, class> class FitnessContainer>
void PathwaysInc<loci>::single_path_finder(
    const intrep::GenericLandscape<loci, F, FitnessContainer>* fitness,
    Path_t path, std::size_t start) {
    path.emplace_back(start);
    bool local_max = true;
    for (std::size_t k = 0; k < loci; k++) {
        if (fitness->at(path.back()) < fitness->at(path.back() ^ (1 << k))) {
            single_path_finder(fitness, path, path.back() ^ (1 << k));
            local_max = false;
        }
    }
    if (local_max) {
        GenericPathways<loci>::pathways.emplace_back(path);
    }
}

template <std::size_t loci>
template <class F, template <std::size_t, class> class FitnessContainer>
void PathwaysInc<loci>::createPathways(
    const intrep::GenericLandscape<loci, F, FitnessContainer>& fitness) {
    GenericPathways<loci>::pathways.clear();
    Path_t path;
    single_path_finder(&fitness, path, 0);
}

template <std::size_t loci>
template <class F, template <std::size_t, class> class FitnessContainer>
void PathwaysInc01<loci>::createPathways(
    const intrep::GenericLandscape<loci, F, FitnessContainer>& fitness) {
    GenericPathways<loci>::pathways.clear();
    Path_t path;
    single_path_finder(&fitness, path, 0);
}

template <std::size_t loci>
template <class T>
void PathwaysInc01<loci>::single_path_finder(T* fitness, Path_t path,
                                             std::size_t start) {
    path.emplace_back(start);
    if (start == (intrep::verticies(loci) - 1)) {
        this->pathways.emplace_back(path);
    } else {
        for (std::size_t k = 0; k < loci; k++) {
            std::size_t neighbour = start | (1 << k);
            if (neighbour != start) {
                if (fitness->at(start) < fitness->at(neighbour)) {
                    single_path_finder(fitness, path, neighbour);
                }
            }
        }
    }
}
} // end namespace intrep
