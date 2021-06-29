#ifndef SLICTOOLS_H
#define SLICTOOLS_H

#include <Imagine/Images.h>

#include "superpixel.h"

bool is_in(const Imagine::Coords<2>&, const Imagine::Image<Imagine::Color>&);

/// The SLIC algorithm
std::vector<Superpixel> SLIC(const Imagine::Image<Imagine::Color>& Img,
                             Imagine::Image<int>& l, int m, int K);

/// Ensure that the Superpixels are connected
void ConnectivityEnforcement(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img);

#endif // SLICTOOLS_H
