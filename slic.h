/**
 * SPDX-License-Identifier: GPL-2.0-or-later
 * @file slic.h
 * @brief SLIC algorithm
 *
 * Copyright (c) 2021 Robin Gay, Pascal Monasse
 * All rights reserved.
 */

#ifndef SLIC_H
#define SLIC_H

#include <Imagine/Images.h>
#include <set>
#include <vector>

/// Are coordinates in the image?
template <typename T>
bool is_in(const Imagine::Coords<2>& p, const Imagine::Image<T>& I)
{ return (0<=p.x() && p.x()<I.width() && 0<=p.y() && p.y() < I.height()); }

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
