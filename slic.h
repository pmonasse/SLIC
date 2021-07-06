/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file slic.h
 * @brief SLIC algorithm
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

#ifndef SLIC_H
#define SLIC_H

#include <Imagine/Images.h>

#include "superpixel.h"
#include <set>
#include <vector>

/// Are coordinates in the image?
template <typename T>
bool is_in(const Imagine::Coords<2>& p, const Imagine::Image<T>& I)
{ return (0<=p.x() && p.x()<I.width() && 0<=p.y() && p.y() < I.height()); }

/// The SLIC algorithm
std::vector<Superpixel> SLIC(const Imagine::Image<Imagine::Color>& Img,
                             Imagine::Image<int>& l, float m, int K);

/// Ensure that the Superpixels are connected
void enforceConnectivity(std::vector<Superpixel>& sp, Imagine::Image<int>& l,
                         const Imagine::Image<Imagine::Color>& Img);

/// Fill set of pixels for each superpixel.
Imagine::Coords<2>* fillSuperpixels(std::vector<Superpixel>& sp,
                                    const Imagine::Image<int>& l);

/// Compute the (compressed) adjacency matrix of superpixels.
void adjacencySuperpixels(const Imagine::Image<int>& l,
                          std::vector< std::set<int> >& adjMatrix);

#endif // SLICTOOLS_H
