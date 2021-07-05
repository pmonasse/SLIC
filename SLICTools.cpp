#include "SLICTools.h"
#include <stack>
#include <queue>
#include <chrono>

bool is_in(int i, int j, const Imagine::Image<Imagine::Color>& img) {
    return(is_in(Imagine::Coords<2>(i,j), img));
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
Imagine::Coords<2> neighbor(int i, int j, int n) {
    assert(n >= 0 && n < 4);
    switch(n) {
        case 0:
            return Imagine::Coords<2>(i+1, j);
        case 1:
            return Imagine::Coords<2>(i, j-1);
        case 2:
            return Imagine::Coords<2>(i-1, j);
        case 3:
            return Imagine::Coords<2>(i, j+1);
        default:
            return Imagine::Coords<2>(i+1, j);
    }
}

/// Return the n-th neighbor of p, starting with right (n=0)
/// and rotating in trigonometric sense
/// \param p Pixel
/// \param n Index of neighbor
/// \return neighbor of p (beware, may be outside image)
Imagine::Coords<2> neighbor(Imagine::Coords<2> p, int n) {
    return neighbor(p.x(), p.y(), n);
}

inline int sq(int x) { return x*x; }

/// Return the squared Euclidian distance between the two colors.
/// \param c1
/// \param c2
/// \return Squared distance
int color_dist(Imagine::Color c1, Imagine::Color c2) {
    return sq(c1.r()-c2.r()) + sq(c1.g()-c2.g()) + sq(c1.b()-c2.b());
}

/// Square norm of gradient at pixel \a p of color image \a Img.
/// The sum of the squares of (horizontal and vertical) derivatives of channels.
int sqNormGradient(const Imagine::Image<Imagine::Color>& Img,
                   const Imagine::Coords<2>& p) {
    int g=0;
    Imagine::Coords<2> q;
    q = p;
    ++q.x();
    if(! is_in(q,Img))
        q.x() = p.x()-1;
    if(! is_in(q,Img))
        q.x() = p.x(); // case w=1
    g += color_dist(Img(q),Img(p));

    q = p;
    ++q.y();
    if(! is_in(q,Img))
        q.y() = p.y()-1;
    if(! is_in(q,Img))
        q.y() = p.y(); // case h=1
    g += color_dist(Img(q),Img(p));
    return g;
}

/// Init superpixels. Their number is adjusted according to the image size.
/// \param[in,out] K Required number of superpixels
/// \param[out] S Dimension of superpixels
/// \param Img Input image
void initSuperpixels(std::vector<Superpixel>& sp, int& K, int& S,
                     const Imagine::Image<Imagine::Color>& Img) {
    const int w=Img.width(), h=Img.height();
    S = (int) sqrt(w*h/(double)K); // Initial size of superpixels

    const int nx = std::max(1,w/S), ny = std::max(1,h/S);
    const int padw = std::max(0,w-S*nx), padh = std::max(0,h-S*ny), s=S/2;
    // Initialization of the superpixels with the color of their center
    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++) {
            int ii=i*S+s+padw/2, jj=j*S+s+padh/2;
            if(is_in(ii,jj,Img))
                sp.push_back(Superpixel(ii, jj, Img(ii,jj), 0));
        }
    K = (int)sp.size();
}

/// Relocalize superpixel centers to lowest gradient position.
/// \param sp Superpixels
/// \param Img Input image
/// \param radius Radius of max motion
/// A negative value of radius is a no-op.
void moveMinGradient(std::vector<Superpixel>& sp,
                     const Imagine::Image<Imagine::Color>& Img, int radius) {
    size_t K = sp.size();
    for(size_t k=0; k<K; k++) {
        int minNorm = std::numeric_limits<int>::max();
        const int x=sp[k].get_x(), y=sp[k].get_y();
        for(int j=-radius; j<=radius; j++)
            for(int i=-radius; i<=radius; i++) {
                Imagine::Coords<2> p(x+i,y+j);
                if(is_in(p, Img)) {
                    int g = sqNormGradient(Img,p);
                    if(g < minNorm) {
                        sp[k].set_x(p.x());
                        sp[k].set_y(p.y());
                        minNorm = g;
                    }
                }
            }
    }
}

/// Assignment step: each pixel is assigned to nearest superpixel
/// \param sp Superpixels
/// \param Img Input image
/// \param m Compactness parameter
/// \param S Influence radius of superpixels
/// \param l Index map of superpixel assignment
/// \param d Distance map in (space,color)
void assignmentStep(const std::vector<Superpixel>& sp,
                    const Imagine::Image<Imagine::Color>& Img,
                    int m, int S, Imagine::Image<int>& l,
                    Imagine::Image<float>& d) {
    int K = (int)sp.size();
    for(int k=0; k<K; k++)
        for(int i=-S; i<S; i++)
            for(int j=-S; j<S; j++) {
                // Look at pixels in a 2S x 2S surrounding of the superpixel
                int ip = sp[k].get_x()+i;
                int jp = sp[k].get_y()+j;

                if(is_in(ip,jp, Img)) {
                    Imagine::Color col = Img(ip,jp);
                    float dist = sp[k].howFar(col, ip,jp, m, S);
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
                const Imagine::Image<int>& l,
                const Imagine::Image<Imagine::Color>& Img) {
    int K = (int)centers.size();
    for(int k=0; k<K; k++)
        for(int l=0; l<6; l++)
            centers[k][l] = 0;

    const int w=Img.width(), h=Img.height();
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            // Adding (x,y,r,g,b) of pixel to its superparent's center
            int pix[5] = {i,j, Img(i,j).r(),Img(i,j).g(),Img(i,j).b()};
            std::vector<int>& c = centers[l(i,j)];
            for(int s=0; s<5; s++)
                c[s] += pix[s];
            ++c[5];
        }

    // Computing barycenters
    for(int k=0; k<K; k++)
        for(int s=0; s<5; s++)
            centers[k][s] /= centers[k][5];
}

/// Squared Euclidian space distance of superpixels' motion
/// \param Superpixels
/// \param centers
float computeError(const std::vector<Superpixel>& sp,
                   const std::vector<std::vector<int>>& centers) {
    float E=0;
    size_t K = sp.size();
    assert(centers.size() == K);
    for(int k=0; k<K; k++)
        E += (centers[k][0]-sp[k].get_x())*(centers[k][0]-sp[k].get_x()) +
             (centers[k][1]-sp[k].get_y())*(centers[k][1]-sp[k].get_y());
    return E;
}

/// Final process of update step: move superpixel centers to their new positions
/// \param sp Superpixels
/// \param centers 6-vector (x,y,r,g,b,n) with n the number of pixels
void moveCenters(std::vector<Superpixel>& sp,
                 const std::vector<std::vector<int>>& centers) {
    size_t K = sp.size();
    for(int k=0; k<K; k++) {
        sp[k].set_x(centers[k][0]);
        sp[k].set_y(centers[k][1]);
        Imagine::Color col(centers[k][2], centers[k][3], centers[k][4]);
        sp[k].set_color(col);
        sp[k].set_size(centers[k][5]);
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
std::vector<Superpixel> SLIC(const Imagine::Image<Imagine::Color>& Img,
                             Imagine::Image<int>& l, int m, int K) {
    auto t0 = std::chrono::high_resolution_clock::now();

    const int w=Img.width(), h=Img.height();

    std::vector<Superpixel> sp;
    int S;  // S will the size of the grid step

    initSuperpixels(sp, K, S, Img);
    l.fill(-1);
    Imagine::Image<float> d(w,h);  // distance map to superpixel's color
    d.fill(float(INFINITY));

    // The minimal gradient initialization doesn't seem to make a difference
    moveMinGradient(sp, Img, 5);

    std::vector<int> il = {0, 0, 0, 0, 0, 0};
    std::vector< std::vector<int> > centers(sp.size(), il);

    float E = 1.0;
    std::cout << "Motions:";
    for(int i=0; E>0; i++) { // Main loop
        assignmentStep(sp, Img, m, S, l, d);
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
void labelCC(Imagine::Image<int>& cc, std::vector<int>& compSizes,
             const Imagine::Image<int>& l) {
    cc.fill(-1);
    std::stack< Imagine::Coords<2> > S;
    const int w=l.width(), h=l.height();
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            if(cc(i,j) != -1)
                continue;
            S.push(Imagine::Coords<2>(i,j));
            int label=l(i,j);
            int labelcc = static_cast<int>(compSizes.size());
            cc(i,j)=labelcc;
            compSizes.push_back(0);
            while(! S.empty()) {
                ++compSizes.back();
                Imagine::Coords<2> p = S.top();
                S.pop();
                for(int n=0; n<4; n++) { // Testing neighbors
                    Imagine::Coords<2> q = neighbor(p,n);
                    if(is_in(q,l) && cc(q)==-1 && l(q)==label) {
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
void discardMinorCC(Imagine::Image<int>& cc,
                    const std::vector<int>& compSizes,
                    Imagine::Image<int>& l, size_t K) {
    const int w=l.width(), h=l.height();
    std::vector<int> maxSizeCC(K,-1); // superpixel -> label of largest cc

    // Fill maxSizeCC
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            int labelcc=cc(i,j);
            int labelS=l(i,j); // Label of superpixel
            int& s = maxSizeCC[labelS];
            if(s<0 || compSizes[s]<compSizes[labelcc])
                s = labelcc;
        }

    // Make orphans the minor cc of each superpixel
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            if(cc(i,j) != maxSizeCC[l(i,j)]) {
                cc(i,j) = -1;
                l(i,j) = -1;
            }
        }
}

/////////////////////////
/// \brief Recomputes the color of the Superpixels as the barycenter of the color of their affiliated pixels
/// /!\ It is assumed that ConnectedComponents already contains only the biggest connected component of each Superpixel,
/// along with orphan pixels (ConnectedComponents(i, j) = -1)
/// \param Superpixels
/// \param K
/// \param ConnectedComponents
/// \param compSizes
/// \param l
/// \param w
/// \param h
/// \param Img
void RecomputeSuperpixelColors(std::vector<Superpixel> Superpixels, int K, Imagine::Image<int> ConnectedComponents, std::vector<int> compSizes, Imagine::Image<int> l, int w, int h, Imagine::Image<Imagine::Color> Img) {

    // newColors will contain the new colors
    Imagine::Image<float> newColors(K, 3);
    for(int k = 0; k < K; k++) {
        for(int s = 0; s < 3; s++) {
            newColors(k, s) = 0;
        }
    }

    // Update newColors by adding the colors of the pixels contained in a Superpixel, divided by the new size of this Superpixel
    // i. e. the size of its biggest connected component
    for(int i = 0; i < w; i ++) {
        for(int j = 0; j < h; j++) {
            if(ConnectedComponents(i, j) != -1) {
                newColors(l(i, j), 0) += float(Img(i, j).r())/float(compSizes[ConnectedComponents(i, j)]);
                newColors(l(i, j), 1) += float(Img(i, j).g())/float(compSizes[ConnectedComponents(i, j)]);
                newColors(l(i, j), 2) += float(Img(i, j).b())/float(compSizes[ConnectedComponents(i, j)]);
            }
        }
    }

    // Update the Superpixels with their new color
    for(int k = 0; k < K; k++) {
        Superpixels[k].set_color(Imagine::Color(int(newColors(k, 0)), int(newColors(k, 1)), int(newColors(k, 2))));
    }
}

/// Assign orphan pixels to nearest (in color) superpixel they are connected to
/// \param sp The superpixels
/// \param[in,out] l
/// \param Img
void assignOrphans(const std::vector<Superpixel>& sp, Imagine::Image<int>& l,
                   const Imagine::Image<Imagine::Color>& Img) {
    const int w=l.width(), h=l.height();

    std::queue< Imagine::Coords<2> > Q;
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
                        Imagine::Coords<2> q = neighbor(i,j,n);
                        if(is_in(q,l) && l(q)<dist) {
                            l(q) = dist;
                            Q.push(q);
                        }
                    }
        if(Qsize == Q.size())
            break;
        Qsize = Q.size();
    }

    while(! Q.empty()) {
        Imagine::Coords<2> p = Q.front();
        Q.pop();
        int nearest = -1;
        int minDist = std::numeric_limits<int>::max();

        for(int n=0; n<4; n++) { // Testing neighbors
            Imagine::Coords<2> q = neighbor(p,n);
            if(! is_in(q,l) || l(q)<0) continue;
            int dist = color_dist(Img(p), sp[l(q)].get_color());
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
void enforceConnectivity(std::vector<Superpixel>& sp, Imagine::Image<int>& l,
                         const Imagine::Image<Imagine::Color>& Img) {
    Imagine::Image<int> cc(l.width(),l.height()); // Labels of cc
    std::vector<int> compSizes; // Sizes of the cc

    auto t0 = std::chrono::high_resolution_clock::now();
    labelCC(cc, compSizes, l);

    discardMinorCC(cc, compSizes, l, sp.size());

    // Recomputing the Superpixel colors taking into account only the biggest connected component (necessary ? NO)
    // RecomputeSuperpixelColors(Superpixels, K, cc, compSizes, l, w, h, Img);

    assignOrphans(sp, l, Img);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time CC: " << std::chrono::duration<double>(t1-t0).count() << std::endl;
}
