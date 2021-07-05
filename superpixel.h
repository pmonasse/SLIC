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

class Superpixel {
public:
    int x, y;                ///< Barycenter coordinates
    int size;                ///< Superpixel size (number of pixels contained)
    Imagine::Color col;      ///< Average color
    Imagine::Coords<2>* pix; ///< Pixels composing the superpixel
public:
    Superpixel(int x0, int y0, const Imagine::Color& c);

    /// \brief Computes the modified R^5 distance (position + color) between this Superpixel's center and the pixel pixel, located
    // in (i, j) in the Imagine::Image<Imagine::Color> img being sliced
    /// \param pixel
    /// \param i
    /// \param j
    /// \param m
    /// \param S
    /// \return float
    float howFar(const Imagine::Color& pixel, int i, int j, int m, int S) const;
    // This distance is the D distance from the article : sqrt((m*euclidianDist/maxEuclidian)^2 + colorDist^2)
    // m is the compactness parameter (user parameter) (weight of the color distance vs. the spatial distance).
    // S is the initial side length of the Superpixel (deduced from the number of Superpixel wished)
};

#endif // SUPERPIXEL_H
