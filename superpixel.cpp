//// ---------------------------------------------------------------------------------------------------------- ////
// PIR - Superpixels (Atelier de Programmation with Pascal Monasse)                                              ///
// superpixel.h                                                                                                  ///
// 22/05/2021 update by Robin GAY                                                                                ///
//                                                                                                               ///
// ----------------------                                                                                        ///
///                                                                                                              ///
///                                                                                                              ///
//// ---------------------------------------------------------------------------------------------------------- ////

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
