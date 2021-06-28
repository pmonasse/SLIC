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

// Superpixel class methods

// Constructors
Superpixel::Superpixel() {}

Superpixel::Superpixel(int x, int y, Imagine::Color col, int sz) {
    // Construct the Superpixel with position, color and size
    _x = x;
    _y = y;
    _col = col;
    _sz = sz;
    _contents = std::vector<Imagine::Coords<2>>(sz);
}

Superpixel::Superpixel(int x, int y, int sz) {
    // Construct a black Superpixel with position and size
    _x = x;
    _y = y;
    _sz = sz;
    _col = Imagine::BLACK;
}

Superpixel::Superpixel(const Superpixel &S) {
    // Copy constructor
    _x = S._x;
    _y = S._y;
    _col = S._col;
    _sz = S._sz;
}

// Destructor
Superpixel::~Superpixel() {}

// Getters
int Superpixel::get_x() const{
    return _x;
}

int Superpixel::get_y() const{
    return _y;
}

Imagine::Color Superpixel::get_color() const{
    return _col;
}

int Superpixel::size() const{
    return _sz;
}

std::vector<Imagine::Coords<2>> Superpixel::get_contents() const{
    return _contents;
}

// Setters
void Superpixel::set_x(const int x) {
    _x = x;
}

void Superpixel::set_y(const int y) {
    _y = y;
}

void Superpixel::set_color(const Imagine::Color col) {
    _col = col;
}

void Superpixel::set_size(const int sz) {
    _sz = sz;
}

void Superpixel::set_contents(std::vector<Imagine::Coords<2>> contents) {
    while(!_contents.empty()) {
        _contents.pop_back();
    }
    for(int i = 0; i < int(contents.size()); i++) {
        _contents.push_back(contents[i]);
    }
    _sz = contents.size();
}

// Other methods
float Superpixel::howFar(const Imagine::Color pixel, int i, int j, int m, int S) const {
    float eucldist = sqrt(float((_x - i)*(_x - i)) + float((_y - j)*(_y - j)));
    float colordist = sqrt(float((_col.r() - pixel.r())*(_col.r() - pixel.r())) + float((_col.g() - pixel.g())*(_col.g() - pixel.g())) + float((_col.b() - pixel.b())*(_col.b() - pixel.b())));
    return(sqrt(float(m*m)*(eucldist*eucldist/float(S*S)) + colordist*colordist));
}

void Superpixel::push_contents(Imagine::Coords<2> pixel) {

    _contents.push_back(pixel);
    _sz += 1;
}

void Superpixel::push_contents(int i, int j) {
    push_contents(Imagine::Coords<2>(i, j));
}


void Superpixel::pop_contents(Imagine::Coords<2> pixel) {
    for(int i = 0; i < int(_contents.size()); i++) {
        if(_contents[i] == pixel) {
            _contents.erase(_contents.begin() + i);
            _sz -= 1;
            break;
        }
    }
}

void Superpixel::pop_contents(int i, int j) {
    pop_contents(Imagine::Coords<2>(i, j));
}
