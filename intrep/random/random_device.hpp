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

#ifndef RANDOM_DEVICE_HPP
#define RANDOM_DEVICE_HPP

#include <fstream>
#include <functional>
#include <chrono>
#include <iostream>

// This is only needed, when std::random_device{} does not work

namespace intrep{
template<class S>
S seed() {
    S random_seed;
    std::ifstream rand_dev("/dev/urandom", std::ios::binary);
    if(rand_dev.is_open()) {
        char* memblock;
        memblock = new char[sizeof(S)];
        rand_dev.read(memblock, sizeof(S));
        rand_dev.close();
        random_seed = *reinterpret_cast<S*>(memblock);
        delete[] memblock;
    }
    else {
        std::cerr << "Chrono" << std::endl;
        auto seed_chrono = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        random_seed = std::hash<S>{}(seed_chrono);
    }
    return random_seed;
}

std::uint64_t seed64() {
    return seed<std::uint64_t>();
}

std::uint32_t seed32() {
    return seed<std::uint32_t>();
}

}
#endif
