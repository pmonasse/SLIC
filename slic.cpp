/**
 * SPDX-License-Identifier: GPL-2.0-or-later
 * @file slic.cpp
 * @brief SLIC algorithm
 *
 * Copyright (c) 2021 Robin Gay, Pascal Monasse
 * All rights reserved.
 */

#include "slic.h"
#include <algorithm>
#include <stack>
#include <queue>
#include <chrono>
#include <iostream>
#include <cmath>
#include <cassert>

inline int sq(int x) { return x*x; }

/// Squared Euclidian distance between the two colors.
/// \param c1 Color
/// \param c2 Color
/// \return Squared distance
int color_dist(const Color& c1, const Color& c2) {
    return sq(c1.r-c2.r) + sq(c1.g-c2.g) + sq(c1.b-c2.b);
}

/// Construct the Superpixel with position, color and size
Superpixel::Superpixel(int x0, int y0, const Color& c) {
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
float Superpixel::dist5D(int i, int j, const Color& c, float wSpace) const {
    int eucldist = sq(x-i)+sq(y-j);
    int colordist = color_dist(col,c);
    return wSpace*wSpace*(float)eucldist + (float)colordist;
}

/// Return the n-th neighbor (4-conn) of (i,j), starting with right (n=0)
/// and rotating in trigonometric sense:
/// i.e. n = 1 -> up neighbor
///      n = 2 -> left neighbor
///      n = 3 -> down neighbor
/// \param i Abscissa
/// \param j Ordinate
/// \param n Index of neighbor
/// \return neighbor of pixel (i,j) (beware, may be outside image)
Pixel neighbor(int i, int j, int n) {
    assert(n >= 0 && n < 4);
    Pixel p = {i,j};
    switch(n) {
    case 0: ++p.x; break;
    case 1: --p.y; break;
    case 2: --p.x; break;
    case 3: ++p.y; break;
    default: assert(false);
    }
    return p;
}

/// Return the n-th neighbor of p, starting with right (n=0)
/// and rotating in trigonometric sense
/// \param p Pixel
/// \param n Index of neighbor
/// \return neighbor of p (beware, may be outside image)
Pixel neighbor(const Pixel& p, int n) {
    return neighbor(p.x, p.y, n);
}

/// Square norm of gradient at pixel \a p of color image \a Img.
/// The sum of the squares of (horizontal and vertical) derivatives of channels.
int sqNormGradient(const Image<Color>& Img, const Pixel& p) {
    int g=0;
    Pixel q = p;
    ++q.x;
    if(! Img.inside(q))
        q.x = p.x-1;
    if(! Img.inside(q))
        q.x = p.x; // case w=1
    g += color_dist(Img(q),Img(p));

    q = p;
    ++q.y;
    if(! Img.inside(q))
        q.y = p.y-1;
    if(! Img.inside(q))
        q.y = p.y; // case h=1
    g += color_dist(Img(q),Img(p));
    return g;
}

/// Init superpixels. Their number is adjusted according to the image size.
/// \param[in,out] K Required number of superpixels
/// \param[out] S Dimension of superpixels
/// \param Img Input image
void initSuperpixels(std::vector<Superpixel>& sp, int& K, int& S,
                     const Image<Color>& Img) {
    const int w=Img.w, h=Img.h;
    S = (int) sqrt(w*h/(double)K); // Initial size of superpixels

    const int nx = std::max(1,w/S), ny = std::max(1,h/S);
    const int padw = std::max(0,w-S*nx), padh = std::max(0,h-S*ny), s=S/2;
    // Initialization of the superpixels with the color of their center
    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++) {
            int ii=i*S+s+padw/2, jj=j*S+s+padh/2;
            if(Img.inside(ii,jj))
                sp.push_back(Superpixel(ii, jj, Img(ii,jj)));
        }
    K = (int)sp.size();
}

/// Relocalize superpixel centers to lowest gradient position.
/// \param sp Superpixels
/// \param Img Input image
/// \param radius Radius of max motion
/// A negative value of radius is a no-op.
void moveMinGradient(std::vector<Superpixel>& sp,
                     const Image<Color>& Img, int radius) {
    size_t K = sp.size();
    for(size_t k=0; k<K; k++) {
        int minNorm = std::numeric_limits<int>::max();
        const int x=sp[k].x, y=sp[k].y;
        for(int j=-radius; j<=radius; j++)
            for(int i=-radius; i<=radius; i++) {
                Pixel p(x+i,y+j);
                if(Img.inside(p)) {
                    int g = sqNormGradient(Img,p);
                    if(g < minNorm) {
                        sp[k].x = p.x;
                        sp[k].y = p.y;
                        minNorm = g;
                    }
                }
            }
    }
}

/// Assignment step: each pixel is assigned to nearest superpixel
/// \param sp Superpixels
/// \param Img Input image
/// \param wSpace Compactness parameter (weight of spatial distance)
/// \param S Influence radius of superpixels
/// \param l Index map of superpixel assignment
/// \param d Distance map in (space,color)
void assignmentStep(const std::vector<Superpixel>& sp,
                    const Image<Color>& Img,
                    float wSpace, int S, Image<int>& l,
                    Image<float>& d) {
    int K = (int)sp.size();
    for(int k=0; k<K; k++) // Look in 2Sx2S surrounding of superpixel center
        for(int i=-S; i<S; i++)
            for(int j=-S; j<S; j++) {
                int ip=sp[k].x+i, jp=sp[k].y+j;
                if(Img.inside(ip,jp)) {
                    Color col = Img(ip,jp);
                    float dist = sp[k].dist5D(ip, jp, col, wSpace);
                    if(d(ip, jp) > dist) {
                        d(ip,jp) = dist;
                        l(ip,jp) = k;
                    }
                }
            }
}

/// First part of update step: returns center of each cluster.
/// \param[out] centers: 6-vector (x,y,r,g,b,n) with n the number of pixels.
/// \param l Index map of superpixels
/// \param Img Input image
void updateStep(std::vector<std::vector<int>>& centers,
                const Image<int>& l, const Image<Color>& Img) {
    int K = (int)centers.size();
    for(int k=0; k<K; k++)
        for(int l=0; l<6; l++)
            centers[k][l] = 0;

    const int w=Img.w, h=Img.h;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            // Adding (x,y,r,g,b) of pixel to its superparent's center
            int pix[5] = {i,j, Img(i,j).r,Img(i,j).g,Img(i,j).b};
            if(l(i,j)<0)
                continue;
            std::vector<int>& c = centers[l(i,j)];
            for(int s=0; s<5; s++)
                c[s] += pix[s];
            ++c[5];
        }

    // Computing barycenters
    for(int k=0; k<K; k++)
        if(centers[k][5]>0)
            for(int s=0; s<5; s++)
                centers[k][s] /= centers[k][5];
}

/// Squared Euclidian space distance of superpixels' motion
/// \param Superpixels
/// \param centers
float computeError(const std::vector<Superpixel>& sp,
                   const std::vector< std::vector<int> >& centers) {
    float E=0;
    size_t K = sp.size();
    assert(centers.size() == K);
    for(int k=0; k<K; k++)
        E += (centers[k][0]-sp[k].x)*(centers[k][0]-sp[k].x) +
             (centers[k][1]-sp[k].y)*(centers[k][1]-sp[k].y);
    return E;
}

/// Final process of update step: move superpixel centers to their new positions
/// \param sp Superpixels
/// \param centers 6-vector (x,y,r,g,b,n) with n the number of pixels
void moveCenters(std::vector<Superpixel>& sp,
                 const std::vector< std::vector<int> >& centers) {
    size_t K = sp.size();
    for(int k=0; k<K; k++) {
        sp[k].x = centers[k][0];
        sp[k].y = centers[k][1];
        Color col(centers[k][2], centers[k][3], centers[k][4]);
        sp[k].col = col;
        sp[k].size = centers[k][5];
    }
}

/// The main SLIC loop.
/// The number of superpixels is deduced from K and the image dimensions.
/// \a l is a map giving for each pixel the index of the assigned superpixel.
/// \param Img The input image
/// \param[out] l Label map of superpixel index
/// \param m Compactness parameter
/// \param K Required number of superpixels
/// \return Collection of superpixels
std::vector<Superpixel> SLIC(const Image<Color>& Img,
                             Image<int>& l, float m, int K) {
    auto t0 = std::chrono::high_resolution_clock::now();

    const int w=Img.w, h=Img.h;

    std::vector<Superpixel> sp;
    int S;  // the size of the grid step

    initSuperpixels(sp, K, S, Img);
    l.fill(-1);
    Image<float> d(w,h);  // distance map to superpixel's color
    d.fill(float(INFINITY));
    float wSpace = m/(float)S; // spatial weight in 5D distance

    moveMinGradient(sp, Img, 5);

    std::vector<int> il = {0, 0, 0, 0, 0, 0};
    std::vector< std::vector<int> > centers(sp.size(), il);

    float E = 1.0;
    std::cout << "Motions:";
    for(int i=0; E>0; i++) { // Main loop
        assignmentStep(sp, Img, wSpace, S, l, d);
        updateStep(centers, l, Img); // Compute new centers of superpixels
        E = computeError(sp, centers);
        moveCenters(sp, centers); // Assign new centers to superpixels
        std::cout << ' ' << E << std::flush;
    }
    std::cout << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t = t1-t0;
    std::cout << "Elapsed time for SLIC: " << t.count() << std::endl;

    return sp;
}

/// Connected components labelling algorithm.
/// \param[out] cc Pixels get label of their connected component.
/// \param[out] compSizes For each label of cc, its number of pixels.
/// \param[in] l Image whose cc of isolevels must be extracted.
void labelCC(Image<int>& cc, std::vector<int>& compSizes, const Image<int>& l) {
    cc.fill(-1);
    std::stack<Pixel> S;
    const int w=l.w, h=l.h;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            if(cc(i,j) != -1)
                continue;
            S.push(Pixel(i,j));
            int label=l(i,j);
            int labelcc = static_cast<int>(compSizes.size());
            cc(i,j)=labelcc;
            compSizes.push_back(0);
            while(! S.empty()) {
                ++compSizes.back();
                Pixel p = S.top();
                S.pop();
                for(int n=0; n<4; n++) { // Testing neighbors
                    Pixel q = neighbor(p,n);
                    if(l.inside(q) && cc(q)==-1 && l(q)==label) {
                        S.push(q);
                        cc(q) = labelcc;
                    }
                }
            }
            // std::cout << compSizes.back() << ' ' << std::flush;
        }
    // std::cout << std::endl;
}

/// Reduce superpixels to their largest connected component
/// \param[in,out] cc connected component labels
/// \param compSizes Size of each cc
/// \param[in,out] l Labels of superpixels
/// \param K Number of superpixels (number of values in l)
void discardMinorCC(Image<int>& cc, const std::vector<int>& compSizes,
                    Image<int>& l, size_t K) {
    const int w=l.w, h=l.h;
    std::vector<int> maxSizeCC(K,-1); // superpixel -> label of largest cc

    // Fill maxSizeCC
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int labelcc=cc(i,j);
            int labelS=l(i,j); // Label of superpixel
            if(labelS>=0) {
                int& s = maxSizeCC[labelS];
                if(s<0 || compSizes[s]<compSizes[labelcc])
                    s = labelcc;
            }
        }

    // Make orphans the minor cc of each superpixel
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            if(l(i,j)>=0 && cc(i,j) != maxSizeCC[l(i,j)]) {
                cc(i,j) = -1;
                l(i,j) = -1;
            }
        }
}

/// Recompute the color of the superpixels.
/// \param sp Superpixels
/// \param l Index map of superpixels
/// \param Img Input image
void computeSuperpixelColors(std::vector<Superpixel>& sp,
                             const Image<int>& l,
                             const Image<Color>& Img) {
    size_t K = sp.size();
    Image<int> col(K,4); // (R,G,B,count)
    col.fill(0);

    const int w=l.w, h=l.h;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int k = l(i,j);
            if(k>=0) {
                Color c = Img(i,j);
                col(k,0) += c.r;
                col(k,1) += c.g;
                col(k,2) += c.b;
                col(k,3)++;
            }
        }

    for(int k=0; k<K; k++)
        if(col(k,3) > 0) {
            Color c(col(k,0)/col(k,3), col(k,1)/col(k,3), col(k,2)/col(k,3));
            sp[k].col = c;
        }
}

/// Assign orphan pixels to nearest (in color) superpixel they are connected to
/// \param sp The superpixels
/// \param[in,out] l
/// \param Img
void assignOrphans(const std::vector<Superpixel>& sp, Image<int>& l,
                   const Image<Color>& Img) {
    const int w=l.w, h=l.h;

    std::queue<Pixel> Q;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++)
            if(l(i,j)<0)
                l(i,j)=std::numeric_limits<int>::min();

    // Do sucessive dilatations until all orphans are queued
    for(int dist=-1; true; --dist) {
        size_t Qsize = Q.size();
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++)
                if(l(i,j) > dist)
                    for(int n=0; n<4; n++) { // Testing neighbors
                        Pixel q = neighbor(i,j,n);
                        if(l.inside(q) && l(q)<dist) {
                            l(q) = dist;
                            Q.push(q);
                        }
                    }
        if(Qsize == Q.size())
            break;
        Qsize = Q.size();
    }

    while(! Q.empty()) {
        Pixel p = Q.front();
        Q.pop();
        int nearest = -1;
        int minDist = std::numeric_limits<int>::max();

        for(int n=0; n<4; n++) { // Testing neighbors
            Pixel q = neighbor(p,n);
            if(! l.inside(q) || l(q)<0) continue;
            int dist = color_dist(Img(p), sp[l(q)].col);
            if(dist < minDist) {
                minDist = dist;
                nearest = l(q);
            }            
        }
        assert(nearest>=0);
        l(p) = nearest;
    }
}

/// Make superpixels connected by keeping only their major cc.
/// Pixels of minor cc get assigned to adjacent superpixel with closest color.
/// \param[in,out] sp The superpixels
/// \param[in,out] l Labels of superpixels
/// \param Img Image
void enforceConnectivity(std::vector<Superpixel>& sp, Image<int>& l,
                         const Image<Color>& Img) {
    Image<int> cc(l.w,l.h); // Labels of cc
    std::vector<int> compSizes; // Sizes of the cc

    auto t0 = std::chrono::high_resolution_clock::now();
    labelCC(cc, compSizes, l);

    discardMinorCC(cc, compSizes, l, sp.size());

    computeSuperpixelColors(sp, l, Img);

    assignOrphans(sp, l, Img);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time CC: " << std::chrono::duration<double>(t1-t0).count() << std::endl;
}

/// Functor to compare pixel coordinates according to their superpixel.
/// Used to sort pixels in \c fillSuperpixels.
struct CompareLabels {
    const Image<int>* l;
    bool operator()(const Pixel& p, const Pixel& q) const {
        return (*l)(p) < (*l)(q);
    }
};

/// Return an allocated array of pixel coordinates. Each superpixel has a
/// pointer \c pix to its first pixel inside this array, its field \c size
/// indicates the number of consecutive pixels.
Pixel* fillSuperpixels(std::vector<Superpixel>& sp, const Image<int>& l){
    const int w=l.w, h=l.h, n=w*h;
    Pixel* c = new Pixel[n];
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++)
            c[j*w+i] = Pixel(i,j);
    CompareLabels cmp;
    cmp.l = &l;
    std::sort(c, c+n, cmp);

    for(int i=0; i<n;) {
        Pixel k=c[i];
        int j=i+1; // Find range of pixels with same superpixel
        while(j<n && !cmp(k,c[j]))
            j++;
        if(l(k)>=0) { // Should always be satisfied if connectivity enforced
            sp[l(k)].size = j-i;
            sp[l(k)].pix = &c[i];
        }
        i=j;
    }
    return c;
}

/// Compute the (compressed) adjacency matrix of superpixels.
/// adjMatrix[k] indicates the index of superpixels adjacent to superpixel k.
/// 4-connectivity is used.
void adjacencySuperpixels(const Image<int>& l,
                          std::vector< std::set<int> >& adjMatrix) {
    const int w=l.w, h=l.h;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++)
            if(l(i,j)>=0)
                for(int n=0; n<2; n++) { // Testing neighbors
                    Pixel q = neighbor(i,j,n);
                    if(l.inside(q) && l(q)>=0) {
                        int k1=l(i,j), k2=l(q);
                        int M = std::max(k1,k2)+1;
                        if(adjMatrix.size() < M) // Make sure indices exist
                            adjMatrix.insert(adjMatrix.end(),M-adjMatrix.size(),
                                             std::set<int>());
                        adjMatrix[k1].insert(k2);
                        adjMatrix[k2].insert(k1);
                    }
                }
}
