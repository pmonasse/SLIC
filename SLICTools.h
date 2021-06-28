#ifndef SLICTOOLS_H
#define SLICTOOLS_H

#include <iostream>
#include <chrono>
// Chrono is used for evaluating complexity

#include <stack>
#include <queue>

#include <Imagine/Images.h>
#include <Imagine/Graphics.h>

#include "superpixel.h"

#include "superpixel.h"

//// ------------------------------------- </> ------------------------------------- ////
////                                  Main functions                                 ////
//// ------------------------------------- </> ------------------------------------- ////

//////////////
/// \brief Loads the image img, and returns it as an Imagine::Image<Imagine::Color>
/// \param img
/// \param w
/// \param h
/// \return
Imagine::Image<Imagine::Color> LoadImage(const char* img, int&w, int&h);

//////////////////
/// \brief Displays the image Img (of dimensions w*h) in the subwindow subwin of the window W
/// \param Img
/// \param W
/// \param subwin
/// \param w
/// \param h
void DisplayImage(Imagine::Image<Imagine::Color> Img, Imagine::Window W, int subwin, int w, int h);

// ----------------------------------------------- CHECKERS AND TESTS (Test version) ----------------------------------------------- //
void GridCheck(int K, int S, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img, Imagine::Window W, int subwin);

void InitStatusCheck(int S, int K, int w, int h, std::vector<Superpixel> Superpixels);

void ConnCompCheck(Imagine::Image<int> ConnectedComponents, int label_num, int w, int h);

void DisplayConnComp(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels);

// Same as previous but displays the colors in the Image
void DisplayColorConnComp(int K, int w, int h, Imagine::Image<int> l, std::vector<Superpixel> Superpixels, Imagine::Image<Imagine::Color> Img);

// Displays the connected components of the Superpixels in random colors
void DisplayRandColorConnComp(int K, int w, int h, Imagine::Image<int> l);

// Displays the centers
void CentersCheck(int K, std::vector<Superpixel> Superpixels);
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

////////////////////////
/// \brief Saves an Image with the name "SLIC-[imageName]"
/// (Yes, it will be used to save the SLICed Image, good call)
/// /!\ Currently, saves the image in the build directory, we'll have to fix that
/// \param SLICImage
/// \param w
/// \param h
/// \param imageName
void SaveSLICImage(Imagine::Image<Imagine::Color> SLICImage, int w, int h, std::string imageName);

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
