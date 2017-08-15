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

#include "generic_landscape.hpp"
#include <random>
#include "../random.hpp"
#include <chrono> // Can be omitted if all functions use device_random.hpp for seeding

namespace intrep{
/**
 * Monotonous landscape. Constructed such that
 * w(\sigma) = w(0) + (1 - w(0)) \frac{|\sigma|}{L}
 */
template<std::size_t loci, class T = double>
class monotonous_landscape : public generic_landscape<loci, T> {
    public:
        /// Constructor set fitness of wildtype and escapetype
        monotonous_landscape(T wildtype, T escapetype);
};

/**
 * Monotonous RMF landscape. This landscape is still single peaked to
 * fulfill the requirements of a single peaked landscape.
 * w(\sigma) = w(0) + (1 - w(0)) \frac{|\sigma|}{L} + random
 * with random being uniform distributed in the intervall [-c, c) s.th.
 * c < \frac{1 - w(0)}{L}. Hence the Var(random) < 1/3 c^2
 */
template<std::size_t loci, class T = double>
class rmf_single_peak : public generic_landscape<loci, T> {
    private:
        T               _variance;
        std::uint32_t   _seed;
        std::mt19937    _engine;
        std::uniform_real_distribution<> _rand;
        T               _wild_type;
        T               _escape_type;
    public:
        /// Constructor with largest possible variance
        rmf_single_peak(T wildtype, T escapetype);
        /// Constructor with given variance
        rmf_single_peak(T wildtype, T escapetype, T variance);
        void create_new();
};

/// House of cards fitness landscape
template<std::size_t loci, class T = double>
class HoC : public generic_landscape<loci, T> {
    protected:
        std::uint32_t   _seed;
        std::mt19937    _engine;
        std::uniform_real_distribution<> _rand;
    public:
        HoC(T min_random, T max_random);
        // Create new random fitness landscape with the same parameters
        void create_new();
};

/**
 * HoC landscape with final type (11...11) being supercritical. All
 * others are subcritical
 */
template<std::size_t loci, class T = double>
class HoC_sub : public HoC<loci, T> {
    private:
        T       _escape_type;
    public:
        // Constructor with random numbers generated in [offset, 1)
        //HoC_sub(T offset, T escapetype, std::uint32_t seed);
        // Default Constructor
        //HoC_sub() : HoC_sub(0, 1, ) {}
        // Constructor with seed generated from entropy pool
        //HoC_sub(T offset, T escapetype) : HoC_sub(offset, escapetype, seed32()) {}
        // Constructor with given offset and given intervall of the random part
        HoC_sub(T min_random, T max_random, T escapetype);
        // Create new landscape
        void create_new();
};

/**
 * Neutral fitness landscape
 */
template<std::size_t loci, class T = double>
class neutral_landscape : public generic_landscape<loci, T> {
    public:
        ///Constructor with value of plateau and escapetype
        neutral_landscape(T plateau, T escapetype);
};

/**
 * Rough Mt Fujii landscape with uniform random part.
 */
template<std::size_t loci, class T>
class RMF_uniform : public generic_landscape<loci, T> {
    private:
        double              c0, c;
        std::mt19937        engine;
        std::uniform_real_distribution<> rand;
    public:
        RMF_uniform(double _c0, double _c, double k) : c0(_c0), c(_c), engine(seed32()), rand(0,k) {
            if((k > 1 - _c0) ||
                    (_c > (1 - _c0 - k) / double(loci - 1)) ||
                    _c0 > 1 || _c < 0) {
                throw std::domain_error("Parameters do not match the requirements\n");
            }
            else {
                create_landscape();
            }
        }
        void create_landscape() {
            for(std::size_t i = 0; i < this->w.size(); i++) {
                this->w[i] = (T)(c * hammingweight(i) + rand(engine) + c0);
            }
        }
        T& back() {
            return this->w.back();
        }
};

/**
 * Rough Mt Fujii landscape. Works only with float, double, long double
 */
template<std::size_t loci, class T = double>
class RMF_landscape :public generic_landscape<loci, T> {
     public:
        template<class S>
        RMF_landscape(double c, double mean, double stddev, S seed) {
            create_landscape(c, mean, stddev, seed);
        }
        RMF_landscape(double c, double mean, double stddev) {
            create_landscape(c, mean, stddev);
        }
        template<class S>
        bool create_landscape(double c, double mean, double stddev, S seed) {
            std::mt19937_64 engine(seed);
            std::normal_distribution<double> dist(mean, stddev);
            for(std::size_t i = 0; i < this->w.size(); i++) {
                this->w[i] = c * hammingweight(i) + dist(engine);
                if(this->w[i] < 0) {
                    this->w[i] = 0;
                }
            }
            return true;
        }
        bool create_landscape(double c, double mean, double stddev) {
            auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            return create_landscape(c, mean, stddev, seed);
        }
};

template<std::size_t loci, class T = double>
class RMF_landscape_glob_max : public generic_landscape<loci, T> {
    private:
        double c;
        std::mt19937_64 engine;
        std::normal_distribution<double> dist;
        std::uint64_t random_seed;
    public:
        template<class S>
        RMF_landscape_glob_max(double _c, double mean, double stddev, double glob_max, S seed) : c(_c), engine(seed), dist(mean, stddev), random_seed(seed){
            this->w.back() = glob_max;
            create_landscape();
        }
        RMF_landscape_glob_max(double c, double mean, double stddev, double glob_max)
            : RMF_landscape_glob_max(c, mean, stddev, glob_max, seed64()) {}
        bool create_landscape() {
            for(std::size_t i = 0; i < this->w.size() - 1; i++) {
                do{
                    this->w[i] = c * hammingweight(i) + dist(engine);
                    if(this->w[i] < 0) {
                        this->w[i] = 0;
                    }
                }while(this->w[i] >= 1);
            }
            return true;
        }
        std::uint64_t getSeed() { return random_seed; }
        void setSeed(std::uint64_t _seed) {
            random_seed = _seed;
            engine.seed(_seed);
        }
};


/******************
 * Implementations
 ******************/


// monotonous_landscape
template<std::size_t loci, class T>
monotonous_landscape<loci, T>::monotonous_landscape(T wildtype, T escapetype) {
    T step = (1 - wildtype) / loci;
    for(std::size_t i = 0; i < this->w.size(); i++){
        this->w[i] = step * hammingweight(i) + wildtype;
    }
    this->w.back() = escapetype;
}

// rmf_single_peak
template<std::size_t loci, class T>
rmf_single_peak<loci, T>::rmf_single_peak(T wildtype, T escapetype, T variance)
    : _variance(variance), _seed(seed32()), _engine(_seed),
        _rand(-std::sqrt(3 * static_cast<double>(_variance)), std::sqrt(3 * static_cast<double>(_variance))),
        _wild_type(wildtype), _escape_type(escapetype) {
    // The test for suitable variance is problematic! Does not work due
    // to rounding error... perhaps?
/*    if(_variance > step * step / 11) {
        std::stringstream except_text;
        except_text << "The variance must be small enough that the fitness landscape is still single peaked.\n";
        except_text << "Variance is " << variance << "\n";
        except_text << "Required maximal variance is " << step * step / 12 << "\n";
        throw std::domain_error(except_text.str());
    }
    */
    create_new();
}

template<std::size_t loci, class T>
rmf_single_peak<loci, T>::rmf_single_peak(T wildtype, T escapetype)
    : rmf_single_peak(wildtype, escapetype, (static_cast<T>(1)-wildtype)*(static_cast<T>(1)-wildtype)/(static_cast<T>(12*loci*loci))) {}

template<std::size_t loci, class T>
void rmf_single_peak<loci, T>::create_new() {
    T step = (static_cast<T>(1) - _wild_type) / static_cast<T>(loci);
    for(std::size_t i = 0; i < this->w.size() - 1; i++) {
        this->w[i] = _wild_type + step * hammingweight(i) + _rand(_engine);
    }
    this->w.back() = _escape_type;
}

// HoC landscape
template<std::size_t loci, class T>
HoC<loci, T>::HoC(T min_random, T max_random) : _seed(seed32()), _engine(_seed), _rand(static_cast<double>(min_random), static_cast<double>(max_random)) {
    create_new();
}

template<std::size_t loci, class T>
void HoC<loci, T>::create_new() {
    for(auto& vert: this->w) {
        vert = _rand(_engine);
    }
}

// HoC_sub landscape
template<std::size_t loci, class T>
HoC_sub<loci, T>::HoC_sub(T min_random, T max_random, T escapetype) :
    HoC<loci, T>::HoC(min_random, max_random), _escape_type(escapetype) {
    this->w.back()      = escapetype;
}

template<std::size_t loci, class T>
void HoC_sub<loci, T>::create_new() {
    HoC<loci, T>::create_new();
    this->w.back() = _escape_type;
}

// Neutral landscape
template<std::size_t loci, class T>
neutral_landscape<loci, T>::neutral_landscape(T plateau, T escapetype) {
    this->w.fill(plateau);
    this->w.back() = escapetype;
}

} // END namespace intrep
