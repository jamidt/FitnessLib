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

#include "pathweight_adaptive/generic_pathweight.hpp"
#include "pathweight_adaptive/01pathweight.hpp"

/**@defgroup adaptivepathweight Adaptive pathway weights
 * @brief Calculation of pathweights of adaptive walks on a fitness landscape
 *
 * This module calculates the pathweights of an adaptive walk on a
 * fitness landscape. That is, the pathweights in the strong selection,
 * weak mutation limit.
 */
