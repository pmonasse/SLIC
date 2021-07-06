/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file superpixel.h
 * @brief Superpixel class
 *
 * Copyright (c) 2021 Robin Gay, Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SUPERPIXEL_H
#define SUPERPIXEL_H

#include <vector>
#include <Imagine/Common.h>

/// Squared Euclidian distance between the two colors.
int color_dist(const Imagine::Color& c1, const Imagine::Color& c2);

/// The superpixel has a center, a set of pixels and a mean color.
class Superpixel {
public:
    int x, y;                ///< Barycenter coordinates
    int size;                ///< Superpixel size (number of pixels contained)
    Imagine::Color col;      ///< Average color
    Imagine::Coords<2>* pix; ///< Pixels composing the superpixel
public:
    Superpixel(int x0, int y0, const Imagine::Color& c);

    /// R^5 distance (position+color) between the center and the pixel.
    float dist5D(int i, int j, const Imagine::Color& c, float wSpace) const;
};

#endif // SUPERPIXEL_H
