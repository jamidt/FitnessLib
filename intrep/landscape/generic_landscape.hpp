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
#include <array>
#include <bitset>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <stdexcept>
#include <string>

#include "../io.hpp"

namespace intrep {
/**@ingroup landscape
 * @brief Generic landscape class
 *
 * This class is rather used for specialization of different fitness
 * landscapes.
 */
template <std::size_t L, class T,
          template <std::size_t, class> class Container = CubeValues>
class GenericLandscape {
  protected:
    Container<L, T> _fitness;

  public:
    using value_type = T;
    // template<class ...Args>
    // GenericLandscape(Args... a) : _fitness{{a...}} {} // Constructor with
    // initialization list

    ///@name Constructor
    ///@{
    /// Default constructor
    GenericLandscape() {}
    /// Copy constructor
    GenericLandscape(const GenericLandscape<L, T, Container>& cp) : _fitness(cp._fitness) {}
    /// Move constructor
    GenericLandscape(GenericLandscape<L, T, Container>&& cp)
        : _fitness(std::move(cp._fitness)) {}
    /// Copyconstructor from container
    GenericLandscape(const Container<L, T>& arr) : _fitness(arr) {}
    /// List initialization
    GenericLandscape(std::initializer_list<T> list) : _fitness(list) {}
    ///@}

    // Operator = returns underlying array with fitness values
    ///@name Assigment operator
    ///@{
    // GenericLandscape<L, T, Container>& operator=(const Container<L, T>&);
    /// Copy assignment
    GenericLandscape<L, T, Container>& operator=(
        const GenericLandscape<L, T, Container>&);
    /// Move assignment operator
    GenericLandscape<L, T, Container>& operator=(
        GenericLandscape<L, T, Container>&&) = default;
    /// List assignment
    GenericLandscape<L, T, Container>& operator=(std::initializer_list<T>);
    ///@}
    ~GenericLandscape<L, T, Container>() = default;

    ///@name Random access
    ///@{
    T& operator[](std::size_t i) noexcept { return _fitness[i]; }
    const T& operator[](std::size_t i) const noexcept { return _fitness[i]; }
    T& at(std::size_t i) { return _fitness.at(i); }
    const T& at(std::size_t i) const { return _fitness.at(i); }
    // T operator()(std::size_t i) const noexcept { return _fitness[i]; }
    // T cat(std::size_t i) const { return _fitness.at(i); }
    ///@}

    T& back() { return _fitness.back(); }
    const T& back() const { return _fitness.back(); }
    T& front() { return _fitness.front(); }
    const T& front() const { return _fitness.front(); }

    void fill(T val) { _fitness.fill(val); }

    ///@name Size
    ///@{
    constexpr std::size_t size() const noexcept { return verticies(L); }
    constexpr std::size_t loci() const noexcept { return L; }
    ///@}

    using iterator = typename Container<L, T>::iterator;
    using const_iterator = typename Container<L, T>::const_iterator;

    ///@name Iterator
    ///@{
    ///@brief Iterator to the underlying container

    /// Returns iterator at the beginning of fitness array
    iterator begin() { return _fitness.begin(); }
    /// Returns iterator at the end of fitness array
    iterator end() { return _fitness.end(); }
    /// Returns constant iterator at beginning of fitness array
    const_iterator cbegin() const { return _fitness.cbegin(); }
    /// Returns constant iterator at end of fitness array
    const_iterator cend() const { return _fitness.cend(); }
    ///@}

    ///@name Output
    ///@{
    /// Print fitness landscape to file in integer format
    void printToFile(const char*, io::loci = io::loci::integer);
    void printToFile(const std::string filename,
                     io::loci p = io::loci::integer) {
        printToFile(filename.c_str(), p);
    }

    ///@}
    ///@name Input
    ///@{
    /// From file in integer format
    void fromFile(std::string filename, io::loci p = io::loci::integer) {
        fromFile(filename.c_str(), p);
    }
    void fromFile(const char*, io::loci = io::loci::integer);
    ///@}

  private:
    /// Print fitness landscape to file in integer format
    void printToFileInt(const char*);

    /// Print fitness landscape to file in binary format
    void printToFileBin(const char*);

    /// Print fitness landscape to file (integer, binary) format
    void printToFileIntBin(const char*);
    /// Print fitness landscape to file (integer, binary) format

    /// From file in integer format
    void fromFileInt(const char*);

    /// From file in binary format
    void fromFileBin(const char*);

};

template<std::size_t L, class T, template<std::size_t, class>class Container>
std::ostream& operator<<(std::ostream& out, const GenericLandscape<L, T, Container>& fitness) {
    for(std::size_t vertex = 0; vertex != fitness.size(); ++vertex) {
        out << vertex << ' ' << fitness[vertex] << '\n';
    }
    return out;
}

//***** GenericLandscape *********
/*template <std::size_t L, class T, template <std::size_t, class> class
Container>
GenericLandscape<L, T, Container>& GenericLandscape<L, T, Container>::operator=(
    const Container<L, T>& arr) {
    _fitness = arr;
    return *this;
}
*/
template <std::size_t L, class T, template <std::size_t, class> class Container>
GenericLandscape<L, T, Container>& GenericLandscape<L, T, Container>::operator=(
    const GenericLandscape<L, T, Container>& landscape) {
    _fitness = landscape._fitness;
    return *this;
}

template <std::size_t L, class T, template <std::size_t, class> class Container>
GenericLandscape<L, T, Container>& GenericLandscape<L, T, Container>::operator=(
    std::initializer_list<T> list) {
    _fitness = list;
    return *this;
}

template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::printToFile(const char* filename,
                                                    io::loci p) {
    switch (p) {
        case io::loci::integer:
            printToFileInt(filename);
            break;
        case io::loci::binary:
            printToFileBin(filename);
            break;
        case io::loci::integer_binary:
            printToFileIntBin(filename);
            break;
        default:
            throw std::invalid_argument("Wrong format selection");
    }
}

// Print landscape
template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::printToFileInt(const char* filename) {
    std::ofstream data(filename);
    for (std::size_t n = 0; n != _fitness.size(); n++) {
        data << n << " " << _fitness[n] << std::endl;
    }
    data.close();
}

// Print landscape
template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::printToFileBin(const char* filename) {
    std::ofstream data(filename);
    for (std::size_t n = 0; n != _fitness.size(); n++) {
        data << std::bitset<L>(n).to_string() << " " << _fitness[n] << std::endl;
    }
    data.close();
}

// Print landscape (integer, binary) format
template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::printToFileIntBin(
    const char* filename) {
    std::ofstream data(filename);
    for (std::size_t n = 0; n != _fitness.size(); n++) {
        data << "(" << n << ", " << std::bitset<L>(n).to_string() << ") "
             << _fitness[n] << std::endl;
    }
    data.close();
}

template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::fromFile(const char* filename,
                                                 io::loci p) {
    switch (p) {
        case io::loci::integer:
            fromFileInt(filename);
            break;
        case io::loci::binary:
            fromFileBin(filename);
            break;
        default:
            throw std::invalid_argument("Wrong argument");
    }
}

// Landscape from file
// Note this reads from the integer representation
// TODO: MUCH MUCH !!!!
template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::fromFileInt(const char* filename) {
    std::ifstream data(filename);
    std::string line;
    while (std::getline(data, line)) {
        std::istringstream lineparsing(line);
        std::string vertex;
        if (std::getline(lineparsing, vertex, ' ')) {
            std::string value;
            std::getline(lineparsing, value);
            this->_fitness[std::stoul(vertex)] = std::stod(value);
        }
    }
    data.close();
}

// Landscape from file
// Note this reads from the binary representation
template <std::size_t L, class T, template <std::size_t, class> class Container>
void GenericLandscape<L, T, Container>::fromFileBin(const char* filename) {
    std::ifstream data(filename);
    std::string line;
    while (std::getline(data, line)) {
        std::string vertex = line.substr(0, L);
        std::string value = line.substr(L + 1, line.size() - (L + 1));
        _fitness.at(std::bitset<L>(vertex).to_ulong()) = std::stod(value);
    }
    data.close();
}
}  // end namespace intrep
