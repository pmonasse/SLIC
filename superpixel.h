//// ---------------------------------------------------------------------------------------------------------- ////
// PIR - Superpixels (Atelier de Programmation with Pascal Monasse)                                              ///
// superpixel.h                                                                                                  ///
// 22/05/2021 update by Robin GAY                                                                                ///
//                                                                                                               ///
// ----------------------                                                                                        ///
/// [Update details]                                                                                             ///                                                                       ///
///                                                                                                              ///
//// ---------------------------------------------------------------------------------------------------------- ////

#ifndef SUPERPIXEL_H
#define SUPERPIXEL_H

#include <vector>
#include <Imagine/Common.h>

/// Squared Euclidian distance between the two colors.
int color_dist(const Imagine::Color& c1, const Imagine::Color& c2);

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
