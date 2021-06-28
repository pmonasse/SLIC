#ifndef SLICTOOLS_H
#define SLICTOOLS_H

#include <Imagine/Images.h>

#include "superpixel.h"

Imagine::Image<Imagine::Color> LoadImage(const char* name);
void SaveImage(const Imagine::Image<Imagine::Color>& img, const char* name);

//////////////////
/// \brief Displays the image Img (of dimensions w*h) in the subwindow subwin of the window W
/// \param Img
/// \param W
/// \param subwin
/// \param w
/// \param h
void DisplayImage(const Imagine::Image<Imagine::Color>& Img, Imagine::Window W, int subwin);

// ----------------------------------------------- CHECKERS AND TESTS (Test version) ----------------------------------------------- //
void InitStatusCheck(int S, int K, int w, int h, std::vector<Superpixel> Superpixels);

void ConnCompCheck(Imagine::Image<int> ConnectedComponents, int label_num, int w, int h);

void DisplayConnComp(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels);

// Same as previous but displays the colors in the Image
void DisplayColorConnComp(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img);

// Displays the connected components of the Superpixels in random colors
void DisplayRandColorConnComp(int K, int w, int h, Imagine::Image<int> l);

// ----------------------------------------------- CHECKERS AND TESTS (Test version) // END ----------------------------------------------- //

//////////////////////////////////
/// \brief Ensures that the Superpixels are connected
/// \param K
/// \param w
/// \param h
/// \param l
/// \param Superpixels
/// \param Img
void ConnectivityEnforcement(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img);

////////////////////////////
/// \brief Fills ImgDestination with a SLICed image, according to user preferences
/// \param superpixels
/// \param borders
/// \param ImgDestination
/// \param Img
/// \param w
/// \param h
/// \param l
/// \param Superpixels
void MakeSLICImage(bool superpixels, bool borders, Imagine::Image<Imagine::Color>& ImgDestination, Imagine::Image<Imagine::Color> Img, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels);

///////////////////////////
/// \brief Gets the compactness param (m), the number of Superpixels desired (K), whether to show the
/// Superpixel borders (displayBorders) and wh. to show the Superpixels (displaySuperpixels)
/// \param m
/// \param K
/// \param displayBorders
/// \param displaySuperpixels
void GetSLICInputs(int& m, int& K, bool& displayBorders, bool& displaySuperpixels);

/////////////////////
/// \brief Initializes the random number generators with current time
void RandInit();

///////////////////
/// \brief The SLIC algorithm
/// \param Img
/// \param m
/// \param K
/// \param l : l(i, j) will be the index of Img(i, j)'s Superparent in returned vector
/// \return A std::vector containing the Superpixels that segment the image
std::vector<Superpixel> SLIC(Imagine::Image<Imagine::Color> Img, Imagine::Image<int> l, int m, int K, int w, int h);

#endif // SLICTOOLS_H
