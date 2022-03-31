/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * @file main.cpp
 * @brief input/output for SLIC algorithm
 *
 * Copyright (c) 2021 Robin Gay, Pascal Monasse
 */

#include "slic.h"
#include "cmdLine.h"
#include "io_jpg.h"
#include "io_png.h"

const Color WHITE(255,255,255);

/// Extensions for JPEG files
const std::string JPG[]={"jpg","JPG","jpeg","JPEG"};
const size_t NJPG = sizeof(JPG)/sizeof(JPG[0]);
/// Extensions for PNG files
const std::string PNG[]={"png","PNG"};
const size_t NPNG = sizeof(PNG)/sizeof(PNG[0]);

/// File extension
std::string ext(const std::string& str) {
    size_t pos = str.rfind('.');
    if(pos==std::string::npos)
        return std::string();
    return str.substr(pos+1, std::string::npos);
}

/// Load image (jpg or png)
unsigned char* load(const char* fname, size_t& nx, size_t& ny) {
    std::string e = ext(fname);
    if(std::find(JPG, JPG+NJPG, e)!=JPG+NJPG)
        return io_jpg_read_u8_rgb(fname, &nx, &ny);
    if(std::find(PNG, PNG+NPNG, e)!=PNG+NPNG)
        return io_png_read_u8_rgb(fname, &nx, &ny);
    std::cerr << "Unkown extension for file " << fname << std::endl;
    return 0;
}

/// Save image (jpg or png)
bool save(const Image<Color>& img, const char* fname) {
    std::string e = ext(fname);
    unsigned char* data = (unsigned char*)img.data;
    if(std::find(JPG, JPG+NJPG, e)!=JPG+NJPG)
        return (io_jpg_write_u8(fname, data, img.w, img.h, 3)==0);
    if(std::find(PNG, PNG+NPNG, e)!=PNG+NPNG)
        return (io_png_write_u8(fname, data, img.w, img.h, 3)==0);
    std::cerr << "Unkown extension for file " << fname << std::endl;
    return false;
}

/// Generate output image
void slic_output(Image<Color>& out,
                 const std::vector<Superpixel>& sp, const Image<int>& l,
                 const Color& col, bool superpixels=true, bool borders=true) {
    const int w=l.w, h=l.h;
    assert(out.w==w && out.h==h);

    // Replacing the color of each pixel by its Superparent's color
    if(superpixels)
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++)
                out(i,j) = l(i,j)>=0? sp[l(i,j)].col: WHITE;

    // Drawing the borders between the superpixels
    if(borders)
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                Pixel n1(i+1,j), n2(i,j+1);
                if((l.inside(n1) && l(n1)!=l(i,j)) ||
                   (l.inside(n2) && l(n2)!=l(i,j)))
                    out(i,j) = col;
            }
}

/// Apply SLIC algorithm and output image.
void slic_image(const Image<Color>& in, Image<Color>& out,
                float m, int K, int g, const Color& col) {
    const int w=in.w, h=in.h;
    Image<int> l(w,h);

    std::vector<Superpixel> sp = SLIC(in, l, m, K, g);
    enforceConnectivity(sp, l, in);

    slic_output(out, sp, l, col);

    // The following calls are useless here, just to show their usage.
    Pixel* pixels = fillSuperpixels(sp, l);
    delete [] pixels;

    std::vector< std::set<int> > adjMatrix;
    adjacencySuperpixels(l, adjMatrix);
}

int main(int argc, char* argv[]) {
    CmdLine cmd;
    int K=1000; // required number of superpixels
    float m=40; // compactness parameter
    Color col = WHITE;
    int g=0;
    cmd.add( make_option('k',K).doc("Required number of superpixels") );
    cmd.add( make_option('m',m).doc("Compactness parameter") );
    cmd.add( make_option('g',g).doc("Radius for minimal gradient search") );
    cmd.add( make_option('c',col).doc("Color of boundary") );
    try { cmd.process(argc, argv);
    } catch(std::string str) {
        std::cerr << "Error: " << str << std::endl;
        return 1;
    }
    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [options] imgIn imgOut\n"
                  << cmd << std::endl;
        return 1;
    }

    size_t nx, ny;
    unsigned char* data = load(argv[1], nx, ny);
    if(! data) {
        std::cerr << "Failed reading image " << argv[1] << std::endl;
        return 1;
    }
    Image<Color> img(nx, ny, data);
    Image<Color> out(img.w, img.h);
    free(data);

    slic_image(img, out, m, K, g, col);

    if(! save(out, argv[2])) {
        std::cerr << "Failed writing image " << argv[2] << std::endl;
        return 1;
    }

    return 0;
}
