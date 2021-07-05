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

// Construct the Superpixel with position, color and size
Superpixel::Superpixel(int x0, int y0, const Imagine::Color& c) {
    x = x0;
    y = y0;
    col = c;
    size = 0;
    pix = 0;
}

// Other methods
float Superpixel::howFar(const Imagine::Color& pixel, int i, int j,
                         int m, int S) const {
    float eucldist = sqrt(float((x - i)*(x - i)) + float((y - j)*(y - j)));
    float colordist = sqrt(float((col.r() - pixel.r())*(col.r() - pixel.r())) + float((col.g() - pixel.g())*(col.g() - pixel.g())) + float((col.b() - pixel.b())*(col.b() - pixel.b())));
    return sqrt(float(m*m)*(eucldist*eucldist/float(S*S))+colordist*colordist);
}
