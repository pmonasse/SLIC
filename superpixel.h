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

class Superpixel
{
    int _x, _y;                                             // Barycenter coordinates
    int _sz;                                                // Superpixel size (number of pixels contained)
    Imagine::Color _col;                                    // Barycenter color
    std::vector<Imagine::Coords<2>> _contents;              // Pixels contained in this Superpixel (necessary ?)
public:
    // Constructors
    Superpixel();
    Superpixel(int x, int y, Imagine::Color col, int sz);   // Additional constructor
    Superpixel(int x, int y, int sz);                       // Other additional constructor
    Superpixel(const Superpixel &S);                        // Copy constructor (uselful ?)

    // Destructor
    ~Superpixel();

    // Getters
    Imagine::Color get_color() const;
    int get_x() const;
    int get_y() const;
    int size() const;
    std::vector<Imagine::Coords<2>> get_contents() const;

    // Setters
    void set_color(const Imagine::Color col);
    void set_x(const int x);
    void set_y(const int y);
    void set_size(const int sz);
    void set_contents(std::vector<Imagine::Coords<2>> contents);

    // Other methods

    ////////////////////
    /// \brief Computes the modified R^5 distance (position + color) between this Superpixel's center and the pixel pixel, located
    // in (i, j) in the Imagine::Image<Imagine::Color> img being sliced
    /// \param pixel
    /// \param i
    /// \param j
    /// \param m
    /// \param S
    /// \return float
    float howFar(Imagine::Color pixel, int i, int j, int m, int S) const;
    // This distance is the D distance from the article : sqrt((m*euclidianDist/maxEuclidian)^2 + colorDist^2)
    // m is the compactness parameter (user parameter) (weight of the color distance vs. the spatial distance).
    // S is the initial side length of the Superpixel (deduced from the number of Superpixel wished)

    ////////////////////
    /// \brief Pushes the pixel pixel in the _contents vector of this Superpixel
    // NB : uses std::vector<Imagine::color>::push_back()
    /// \param pixel
    void push_contents(Imagine::Coords<2> pixel);
    void push_contents(int i, int j);

    //////////////////
    /// \brief Pops the pixel from the Superpixel's contents (no use)
    /// \param pixel
    void pop_contents(Imagine::Coords<2> pixel);
    void pop_contents(int i, int j);
};

#endif // SUPERPIXEL_H
