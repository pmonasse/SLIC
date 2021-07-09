/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * @file image.h
 * @brief Elementary image class
 *
 * Copyright (c) 2021 Robin Gay, Pascal Monasse
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <cstring>

/// Pixel location
struct Pixel {
    int x,y;
    Pixel() {}
    Pixel(int x0, int y0): x(x0), y(y0) {}
};

/// RGB Color
struct Color {
    unsigned char r,g,b;
    Color() {}
    Color(unsigned char r0, unsigned char g0, unsigned char b0)
    : r(r0), g(g0), b(b0) {}    
};

/// A simple (simplistic?) image class
template <typename T>
struct Image {
    const int w,h;
    T* data;
public:
    Image(int w0, int h0);
    Image(int w0, int h0, unsigned char* rawdata);
    ~Image();

    void fill(T value);

    T operator()(int i, int j) const;
    T operator()(const Pixel& p) const {
        return operator()(p.x,p.y);
    }
    T& operator()(int i, int j);
    T& operator()(const Pixel& p) {
        return operator()(p.x,p.y);
    }

    bool inside(int x, int y) const {
        return (0<=x && x<w && 0<=y && y<h);
    }
    bool inside(const Pixel& p) const {
        return inside(p.x, p.y);
    }
private:
    Image(const Image&);
    void operator=(const Image&);
};

/// Constructor
template <typename T>
Image<T>::Image(int w0, int h0)
: w(w0), h(h0) {
    data = new T[w*h];
}

/// Constructor copying data
template <typename T>
Image<T>::Image(int w0, int h0, unsigned char* rawdata)
: w(w0), h(h0) {
    data = new T[w*h];
    std::memcpy(data, rawdata, w*h*sizeof(T));
}

/// Destructor
template <typename T>
Image<T>::~Image() {
    delete [] data;
}

/// Fill with constant value
template <typename T>
void Image<T>::fill(T value) {
    T* p=data;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++)
            *p++ = value;
}

/// Pixel access (read-only)
template <typename T>
T Image<T>::operator()(int i, int j) const {
    return data[i+j*w];
}

/// Pixel access (read/write)
template <typename T>
T& Image<T>::operator()(int i, int j) {
    return data[i+j*w];
}

#endif
