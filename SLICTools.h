#ifndef SLICTOOLS_H
#define SLICTOOLS_H

#include <Imagine/Images.h>

#include "superpixel.h"
#include <set>

/// Are coordinates in the image?
template <typename T>
bool is_in(const Imagine::Coords<2>& p, const Imagine::Image<T>& I)
{ return (0<=p.x() && p.x()<I.width() && 0<=p.y() && p.y() < I.height()); }

/// The SLIC algorithm
std::vector<Superpixel> SLIC(const Imagine::Image<Imagine::Color>& Img,
                             Imagine::Image<int>& l, int m, int K);

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
