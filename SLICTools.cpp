#include "SLICTools.h"
#include <stack>
#include <queue>
#include <chrono>

//// ------------------------------------- </> ------------------------------------- ////
////                                 Useful functions                                ////
//// ------------------------------------- </> ------------------------------------- ////

////////////////
/// \brief If is_in == true, coordinates of pix are in the Image img
/// \param pix
/// \param img
/// \return in
bool is_in(Imagine::Coords<2> pix, Imagine::Image<Imagine::Color> img) {
    return(pix.x() >= 0 && pix.x() < img.width() && pix.y() >= 0 && pix.y() < img.height());
}

////////////////
/// \brief If is_in == true, pixel (i, j) is in img
/// \param i
/// \param j
/// \param img
/// \return in
bool is_in(int i, int j, Imagine::Image<Imagine::Color> img) {
    return(is_in(Imagine::Coords<2>(i, j), img));
}

////////////////////
/// \brief Returns the n-th neighbor (4-conn) of (i, j), starting with right (n = 0) and rotating in trigonometric sense
/// i.e. n = 1 -> up neighbor
///      n = 2 -> left neighbor
///      n = 3 -> down neighbor
/// \param i
/// \param j
/// \param n
/// \return Imagine::Coords<2> neighbor
Imagine::Coords<2> neighbor(int i, int j, int n) {
    assert(n >= 0 && n < 4);
    switch(n) {
        case 0 :
            return Imagine::Coords<2>(i + 1, j);
        case 1 :
            return Imagine::Coords<2>(i, j - 1);
        case 2 :
            return Imagine::Coords<2>(i - 1, j);
        case 3 :
            return Imagine::Coords<2>(i, j + 1);
        default :
            return Imagine::Coords<2>(i + 1, j);
    }
}

////////////////////
/// \brief Returns the n-th neighbor of pix, starting with right (n = 0) and rotating in trigonometric sense
/// \param pix
/// \param n
/// \return Imagine::Coords<2> neighbor
Imagine::Coords<2> neighbor(Imagine::Coords<2> pix, int n) {
    return neighbor(pix.x(), pix.y(), n);
}

/////////////////
/// \brief Returns the squared euclidian distance between the two colors (R, G, B)
/// \param col1
/// \param col2
/// \return int
int color_dist(Imagine::Color col1, Imagine::Color col2) {
    return((col1.r() - col2.r())*(col1.r() - col2.r()) + (col1.g() - col2.g())*(col1.g() - col2.g()) + (col1.b() - col2.b())*(col1.b() - col2.b()));
}

//// ------------------------------------- </> ------------------------------------- ////
////                                  Main functions                                 ////
//// ------------------------------------- </> ------------------------------------- ////

/////////////////
/// \brief Initializes l and d (cf. article), recalculates K to adapt it to the image Img, computes the distance S
/// \param l
/// \param d
/// \param K
/// \param S
/// \param w
/// \param h
/// \param Img
/// \return
void InitGrid(std::vector<Superpixel>& Superpixels,
              Imagine::Image<int>& l, Imagine::Image<float>& d,
              int& K, int& S,
              const Imagine::Image<Imagine::Color>& Img) {
    const int w=Img.width(), h=Img.height();
    S = int(sqrt((w*h)/K));
    // The initial size of the Superpixels (cf. article)

    // Update of K so that the number of Superpixels stored is the actual number, not the one required by the user
    // (because bc of the image dimensions, the exact number required is almost never reached with a regular grid)
    int newK = (w/S) * (h/S);
    K = newK;

    // Initialization of d and l (cf. article)
    for(int i=0; i<w; i++) {
        for(int j=0; j<h; j++) {
            l(i,j) = -1;
            d(i,j) = float(INFINITY);
        }
    }

    // These two int are the difference bewteen the width (resp. height) of the Image and the width (resp. height) of the grid
    // They will help center the grid to avoid missing pixels
    const int padw = w - S*(w/S), padh = h - S*(h/S), s=S/2;
    // Initialization of the superpixels with the color of their center
    for(int j=0; j*S<h; j++)
        for(int i=0; i*S<w; i++)
            Superpixels.push_back(Superpixel(i*S+s+padw/2, j*S+s+padh/2,
                                             Img(i*S+s+padw/2,j*S+s+padh/2),0));
}

///////////////////
/// \brief Does the assignment step according to the article
/// \param sp
/// \param Img
/// \param K
/// \param m
/// \param S
/// \param l
/// \param d
void AssignmentStep(const std::vector<Superpixel>& sp,
                    const Imagine::Image<Imagine::Color>& Img,
                    int K, int m, int S, Imagine::Image<int>& l,
                    Imagine::Image<float>& d) {
    for(int k=0; k<K; k++)
        for(int i=-S; i<S; i++)
            for(int j=-S; j<S; j++) {
                // Looking at the pixels in a 2S x 2S surrounding of the Superpixel sp[k]
                int ip = sp[k].get_x()+i;
                int jp = sp[k].get_y()+j;

                if(is_in(ip, jp, Img)) {
                    Imagine::Color col = Img(ip,jp);
                    float dist = sp[k].howFar(col, ip,jp, m, S);
                    if(d(ip, jp) > dist) {
                        d(ip,jp) = dist;
                        l(ip,jp) = k;
                    }
                }
            }
}

///////////////////
/// \brief Does the first part of the update step and returns a K-element vector of 5-element vectors of int, containing for each cluster center
/// its new color and position coordinates
/// \param K
/// \param w
/// \param h
/// \param l
/// \param Img
/// \return
void UpdateStep(std::vector<std::vector<int>>& centers, int K,
                const Imagine::Image<int>& l,
                const Imagine::Image<Imagine::Color>& Img) {

    for(int k=0; k<K; k++)
        for(int l=0; l<6; l++)
            centers[k][l] = 0;

    // A vector containing in this order the future x, y, r, g, b of the new barycenters
    // and a counter of the number of pixels in the updated Superpixel represented by centers[k]

    const int w=Img.width(), h=Img.height();
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            // Adding the values (x, y, r, g, b) of each pixel to its Superparent center's vector in centers
            int currPix[5] = {i,j, Img(i,j).r(),Img(i,j).g(),Img(i,j).b()};
            for(int s=0; s<5; s++)
                centers[l(i,j)][s] += currPix[s];
            ++centers[l(i,j)][5];
        }

    // Computing the proper barycenter of the pixels of each Superpixels
    // by dividing each value by the number of pixels in the Superpixel
    for(int k=0; k<K; k++)
        for(int s=0; s<5; s++)
            centers[k][s] /= centers[k][5];
}

///////////////////
/// \brief Computes the squared Euclidian distance between the new positions of the superpixels (stocked in centers) and their former positions
/// \param E
/// \param K
/// \param Superpixels
/// \param centers
float ComputeError(int K, const std::vector<Superpixel>& sp,
                  const std::vector<std::vector<int>>& centers) {
    float E=0;
    for(int k=0; k<K; k++)
        // Computing the squared euclidian distance between the previous positions of the centers and the new positions of the centers
        // Maybe we could improve the algorithm wy choosing another distance for error computation ?
        E += (centers[k][0]-sp[k].get_x())*(centers[k][0]-sp[k].get_x()) +
             (centers[k][1]-sp[k].get_y())*(centers[k][1]-sp[k].get_y());
    return E;
}

////////////////////////
/// \brief Does the final process of the update step, i. e. moving the Superpixel centers to their new positions, still stocked in centers
/// \param Superpixels
/// \param centers
/// \param K
void MoveCenters(std::vector<Superpixel>& Superpixels,
                 const std::vector<std::vector<int>>& centers, int K) {
    // Moving the barycenters to their new positions
    for(int k=0; k<K; k++) {
        // Updating the values of the Superpixels
        Superpixels[k].set_x(centers[k][0]);
        Superpixels[k].set_y(centers[k][1]);
        Imagine::Color col(centers[k][2], centers[k][3], centers[k][4]);
        Superpixels[k].set_color(col);
        Superpixels[k].set_size(centers[k][5]);
    }
}

/////////////////////////
/// \brief A basic depth-search connected components algorithm.
/// NB: Connected components in the sense of the affiliation to Superpixels, given by l \n
/// NB: max_label is redirected to the "maximum of the labels", which is ///THE NUMBER OF LABELS +1///
/// \param ConnectedComponents
/// \param max_label
/// \param l
/// \param h
/// \param w
/// \param Img
void ConnectedComponentsAlgorithm(Imagine::Image<int> ConnectedComponents, int& max_label, std::vector<int>& compSizes, Imagine::Image<int> l, int h, int w, Imagine::Image<Imagine::Color> Img) {

    std::stack<Imagine::Coords<2>> pile;
    // The stack used in the connected components algorithm

    int label = 0;

    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h; j++) {

            if(ConnectedComponents(i, j) == -1) {
                ConnectedComponents(i, j) = label;
                pile.push(Imagine::Coords<2>(i, j));

                compSizes.push_back(1);

                while(!pile.empty()) {

//                    // Status check (test version)
//                    std::cout << "Pile size : " << pile.size() << std::endl;

                    Imagine::Coords<2> pix = pile.top();
                    pile.pop();

                    // Testing if the 4 neighbours of Img(pix) are in Img(pix)'s connected component
                    for(int n = 0; n < 4; n++) {
                        Imagine::Coords<2> nei = neighbor(pix, n);
                        if(is_in(nei, Img)) {
                            if(l(nei) == l(i, j) && ConnectedComponents(nei) == -1) {
                                ConnectedComponents(nei) = ConnectedComponents(pix);
                                pile.push(nei);

                                compSizes[ConnectedComponents(i, j)] += 1;
                            }
                        }
                    }
                }
                label += 1;
            }
        }
    }
    max_label = label;
}

////////////////////////////
/// \brief Reduces the Superpixels to their largest connected component and discards all other connected components
/// \param compSizes
/// \param ConnectedComponents
/// \param l
/// \param K
/// \param w
/// \param h
void KeepOnlyMaxSizesComponents(std::vector<int> compSizes, Imagine::Image<int> ConnectedComponents, Imagine::Image<int> l, int K, int w, int h) {

    // maxSizeComponents will contain for each Superpixel k the label of its largest connected subset
    std::vector<int> maxSizeComponents(K);

    // 1st step: determine the label of the largest connected component for each Superpixel k
    // i. e. filling maxSizeComponents
    for(int k = 0; k < K; k++) {
        int max_size = -1;
        int max_label = -1;
        for(int i = 0; i < w; i++) {
            for(int j = 0; j < h; j++) {
                if(l(i, j) == k && compSizes[ConnectedComponents(i, j)] > max_size) {
                    max_size = compSizes[ConnectedComponents(i, j)];
                    max_label = ConnectedComponents(i, j);
                }
            }
        }
        maxSizeComponents[k] = max_label;
    }

    // 2nd step: discard the connected components of the remaining pixels and de-assign them from their Superparent
    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h; j++) {
            if(ConnectedComponents(i, j) != maxSizeComponents[l(i, j)]) {
                ConnectedComponents(i, j) = -1;
                l(i, j) = -1;
            }
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

    // newColors isn't useful anymore
    newColors.~Image();
}

///////////////////////////
/// \brief Reassigns the orphan pixels to the nearest (in color) Superpixel they are connected to
/// \param Superpixels
/// \param l
/// \param w
/// \param h
/// \param Img
void ReassignPixels(std::vector<Superpixel> Superpixels, Imagine::Image<int> l, int w, int h, Imagine::Image<Imagine::Color> Img) {

    // This queue will contain all the orphan pixels in the Image.
    // Using a queue allows to memorize which pixels are isolated and thus cannot be immediately reassigned to a Superpixel they are connected to
    std::queue<Imagine::Coords<2>> queue;

    // Initialization of the queue
    for(int i = 0; i < w; i++) {
        for(int j =0; j < h; j++) {
            if(l(i, j) == -1) {
                queue.push(Imagine::Coords<2>(i, j));
            }
        }
    }

    // While the queue is not empty, it means that there still are isolated pixels that have not been treated yet
    while(!queue.empty()) {

        // Going to treat the front pixel of the queue
        Imagine::Coords<2> pix = queue.front();
        // Pop the front element
        queue.pop();
        // A boolean checking if it is isolated
        bool isolated = true;
        // Look at pix's neighbors

        //// Determine the color-nearest Superpixel that is connected to pix

        // Initializing color distance and closest Superpixel index
        int color_nearest_superpix = -1;
        int min_color_dist = int(INFINITY);

        for(int n = 0; n < 4; n++) {

            Imagine::Coords<2> nei = neighbor(pix, n);

            // Check if the neighbor is in the Image
            if(is_in(nei, Img)) {

                // If l(nei) = -1, it means that nei is not in a Superpixel's connected component
                // Thus, nei is not taken into account
                if(l(nei) == -1) {
                    continue;
                }

                // Else, determine if nei's Superparent is closer in color to pix than the former neighbor treated
                // (or than infinity)
                else {

                    // If there is a neighbor belonging to a Superpixel, it means pix is not isolated
                    isolated = false;

                    // currCD: current color distance
                    int currCD = color_dist(Img(pix), Superpixels[l(nei)].get_color());

                    // If nei's Superpixel is closer to pix in color than the former neighbor,
                    // replace the nearest Superpixel with it
                    if(currCD < min_color_dist) {
                        min_color_dist = currCD;
                        color_nearest_superpix = l(nei);
                    }
                }
            }
        }

        // If at the end of the treatment of the four neighbors of pix, isolated is still true,
        // it means that pix is isolated. It is thus put back in the queue to be treated later
        if(isolated) {
            queue.push(pix);
        }

        // Else, reassign pix to the color-closest Superpixel connected to it
        else {
            l(pix) = color_nearest_superpix;
        }
    }
    // Once the queue is empty, all pixel have been reassigned
}

void ConnectivityEnforcement(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img) {

    Imagine::Image<int> ConnectedComponents(w, h);
    // This Image will contain for each pixel, which connected component it is in
    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h; j++) {
            ConnectedComponents(i, j) = -1;
        }
    }

    // This int will be the number of labels (+1 since numbering starts at 0)
    int label_num;
    // compSizes will be filled with the sizes of the components by the ConnectedComponentsAlgorithm function
    std::vector<int> compSizes;

    ConnectedComponentsAlgorithm(ConnectedComponents, label_num, compSizes, l, h, w, Img);

    KeepOnlyMaxSizesComponents(compSizes, ConnectedComponents, l, K, w, h);

//        // Recomputing the Superpixel colors taking into account only the biggest connected component (necessary ? NO)
//        RecomputeSuperpixelColors(Superpixels, K, ConnectedComponents, compSizes, l, w, h, Img);

    // Reassigning the orphan pixels to the nearest connected Superpixel
    ReassignPixels(Superpixels, l, w, h, Img);
}

//////////////////
/// \brief Computes the euclidian norm of a 2-element FVector of Imagine::Color
/// \param Vect
/// \return
int basicNorm(Imagine::FVector<Imagine::Color, 2> Vect) {
    int a = Vect[0].r()*Vect[0].r() + Vect[0].g()*Vect[0].g() + Vect[0].b()*Vect[0].b();
    int b = Vect[1].r()*Vect[1].r() + Vect[1].g()*Vect[1].g() + Vect[1].b()*Vect[1].b();
    return(a + b);
}

/////////////////////////////
/// \brief Relocalizes the Superpixel centers in the lowest gradient position in a square of zoneSize by zoneSize pixels around their initial positions
/// \param K
/// \param Superpixels
/// \param Img
void InitMinGradient(int K, std::vector<Superpixel>& Superpixels, Imagine::Image<Imagine::Color> Img, int zoneSize) {

    for(int k = 0; k < K; k++) {

        // kGrad is the gradient of the current Superpixel center's position
        Imagine::FVector<Imagine::Color, 2> kGrad = Imagine::gradient(Img, Imagine::Coords<2>(Superpixels[k].get_x(), Superpixels[k].get_y()));

        int minNorm = basicNorm(kGrad);
        int currX = Superpixels[k].get_x();
        int currY = Superpixels[k].get_y();

        // For every pixel in a zoneSize x zoneSize surrounding
        for(int i = -(zoneSize/2); i < zoneSize/2 + 1; i++) {
            for(int j = -(zoneSize/2); j < zoneSize/2 + 1; j++) {

                Imagine::Coords<2> currPix(currX + i, currY + j);

                // Check if the neighbour is in the Image
                if(is_in(currPix, Img)) {

                    // currGrad is the gradient in the current neighbor
                    Imagine::FVector<Imagine::Color, 2> currGrad = gradient(Img, currPix);

                    // If the norm of the gradient of this neighbor is lowest than the lowest so far, replace the lowest by it and move the Superpixel center
                    if(basicNorm(currGrad) < minNorm) {
                        Superpixels[k].set_x(currPix.x());
                        Superpixels[k].set_y(currPix.y());
                        minNorm = basicNorm(currGrad);
                    }
                }
            }
        }
    }
}

void MakeSLICImage(bool superpixels, bool borders,
                   Imagine::Image<Imagine::Color>& ImgDestination,
                   const Imagine::Image<Imagine::Color>& Img,
                   const Imagine::Image<int>& l,
                   const std::vector<Superpixel>& Superpixels) {
    const int w=Img.width(), h=Img.height();
    assert(ImgDestination.width() == w && ImgDestination.height() == h);

    // Replacing the color of each pixel by its Superparent's color
    if(superpixels)
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++)
                ImgDestination(i,j) = Superpixels[l(i,j)].get_color();

    // Drawing the borders between the Superpixels
    if(borders)
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++)
                for(int n=0; n<2; n++) {
                    Imagine::Coords<2> nei = neighbor(i, j, n);
                    if(is_in(nei, Img) && l(nei) != l(i, j))
                        ImgDestination(i, j) = Imagine::WHITE;
                }
}

void GetSLICInputs(int& m, int& K, bool& displayBorders, bool& displaySuperpixels) {
    // Inputs
    std::cout << "Number of Superpixels: ";
    std::cin >> K;
    std::cout << "Compactness parameter (>=1): ";
    std::cin >> m;
    std::cout << "Display borders? ";
    std::cin >> displayBorders;
    std::cout << "Display Superpixels? ";
    std::cin >> displaySuperpixels;
}

std::vector<Superpixel> SLIC(const Imagine::Image<Imagine::Color>& Img,
                             Imagine::Image<int>& l, int m, int K) {

    // Chronometer for test version
    auto time_before_SLICLoops = std::chrono::high_resolution_clock::now();

    // Status check
    std::cout << "Computing SLIC..." << std::endl;
    const int w=Img.width(), h=Img.height();

    std::vector<Superpixel> Superpixels;    // A vector that will contain all the Superpixels
    Imagine::Image<float> d(w,h);  // d(i, j) will be the distance from Img(i, j) to its Superparent Superpixels[l(i, j)]
    int S;  // S will the size of the grid step

    // InitSuperpixels initializes l, d, S and Superpixels
    InitGrid(Superpixels, l, d, K, S, Img);

//    // The minimal gradient initialization doesn't seem to make a difference
//    InitMinGradient(K, Superpixels, Img, 10);

    // centers will be useful for the update step, as it will contain the new positions of the Superpixels,
    // allowing computation of the error (see UpdateStep)
    // /!\ centers MUST be initialized after the grid because K is modified by InitGrid
    std::vector<int> il = {0, 0, 0, 0, 0, 0};
    std::vector<std::vector<int>> centers(K, il);

    float E = 0.;
    int loopCounter = 0;

    // Main loop
    do {
        loopCounter += 1;
        AssignmentStep(Superpixels, Img, K, m, S, l, d);
        UpdateStep(centers, K, l, Img);
        // Now centers contains the position and color coordinates of the new barycenters, and the sizes of the corresponding Superpixel

        E = ComputeError(K, Superpixels, centers);
        MoveCenters(Superpixels, centers, K);
//        std::cout << "Error : " << E << " at loop " << loopCounter << std::endl;;
    } while(E > 0);
    // The stop condition is E <= 0 : since E is an integer and converges towards zero, this condition always makes the algorithm end

    // Chronometer for test version
    auto time_after_SLICLoops = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time_SLICLoops = (time_after_SLICLoops - time_before_SLICLoops);
    std::cout << "Elapsed time for SLIC init and loops: " << elapsed_time_SLICLoops.count() << std::endl;

    return Superpixels;
}
