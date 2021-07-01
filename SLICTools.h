#ifndef SLICTOOLS_H
#define SLICTOOLS_H

#include <Imagine/Images.h>

#include "superpixel.h"

/// Are coordinates in the image?
template <typename T>
bool is_in(const Imagine::Coords<2>& p, const Imagine::Image<T>& I)
{ return (0<=p.x() && p.x()<I.width() && 0<=p.y() && p.y() < I.height()); }

/// The SLIC algorithm
std::vector<Superpixel> SLIC(const Imagine::Image<Imagine::Color>& Img,
                             Imagine::Image<int>& l, int m, int K);

/// Ensure that the Superpixels are connected
void ConnectivityEnforcement(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img);

#endif // SLICTOOLS_H
