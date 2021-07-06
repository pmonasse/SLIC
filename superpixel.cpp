/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file superpixel.cpp
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

#include "superpixel.h"

inline int sq(int x) { return x*x; }

/// Squared Euclidian distance between the two colors.
/// \param c1 Color
/// \param c2 Color
/// \return Squared distance
int color_dist(const Imagine::Color& c1, const Imagine::Color& c2) {
    return sq(c1.r()-c2.r()) + sq(c1.g()-c2.g()) + sq(c1.b()-c2.b());
}

// Construct the Superpixel with position, color and size
Superpixel::Superpixel(int x0, int y0, const Imagine::Color& c) {
    x = x0;
    y = y0;
    col = c;
    size = 0;
    pix = 0;
}

/// R^5 distance (position+color) between the center and the pixel.
/// \param i Abscissa of pixel
/// \param j Ordinate of pixel
/// \param c Color of pixel
/// \param wSpace Compactness parameter (weight of spatial distance)
/// \return Squared distance
float Superpixel::dist5D(int i, int j, const Imagine::Color& c,
                         float wSpace) const {
    int eucldist = sq(x-i)+sq(y-j);
    int colordist = color_dist(col,c);
    return wSpace*wSpace*(float)eucldist + (float)colordist;
}
