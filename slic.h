/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * @file slic.h
 * @brief SLIC algorithm (superpixels)
 *
 * Copyright (c) 2021-2022 Robin Gay, Pascal Monasse
 */

#ifndef SLIC_H
#define SLIC_H

#include "image.h"
#include <set>
#include <vector>

/// The superpixel has a center, a set of pixels and a mean color.
class Superpixel {
public:
    float x,y;  ///< Barycenter coordinates
    int size;   ///< Superpixel size (number of pixels contained)
    Color col;  ///< Average color
    Pixel* pix; ///< Pixels composing the superpixel
public:
    Superpixel(int x0, int y0, const Color& c);

    /// R^5 distance (position+color) between the center and the pixel.
    float dist5D(int i, int j, const Color& c, float wSpace) const;
};

/// The SLIC algorithm
std::vector<Superpixel> SLIC(const Image<Color>& I,
                             Image<int>& l, float m, int K, int g=0);

/// Ensure that the Superpixels are connected
void enforceConnectivity(std::vector<Superpixel>& sp, Image<int>& l,
                         const Image<Color>& I);

/// Fill set of pixels for each superpixel.
Pixel* fillSuperpixels(std::vector<Superpixel>& sp, const Image<int>& l);

/// Compute the (compressed) adjacency matrix of superpixels.
void adjacencySuperpixels(const Image<int>& l,
                          std::vector< std::set<int> >& adjMatrix);

#endif // SLICTOOLS_H
