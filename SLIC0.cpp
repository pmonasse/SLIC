//// ---------------------------------------------------------------------------------------------------------- ////
// PIR - Superpixels (Atelier de Programmation with Pascal Monasse)                                              ///
// 22/05/2021 update by Robin GAY                                                                                ///
//                                                                                                               ///
// ----------------------                                                                                        ///
/// [Update details]                                                                                             ///
/// - First working version resulting into ACTUAL (connex) Superpixels                                           ///
///                                                                                                              ///
//// ---------------------------------------------------------------------------------------------------------- ////

#include "SLICTools.h"
#include <iostream>

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

void MakeSLICImage(bool superpixels, bool borders,
                   Imagine::Image<Imagine::Color>& ImgDestination,
                   const std::vector<Superpixel>& sp,
                   const Imagine::Image<int>& l) {
    const int w=l.width(), h=l.height();
    assert(ImgDestination.width()==w && ImgDestination.height()==h);

    // Replacing the color of each pixel by its Superparent's color
    if(superpixels)
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++)
                ImgDestination(i,j) = l(i,j)>=0? sp[l(i,j)].col: Imagine::WHITE;

    // Drawing the borders between the superpixels
    if(borders)
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                Imagine::Coords<2> n1(i+1,j), n2(i,j+1);
                if((is_in(n1,l) && l(n1)!=l(i,j)) ||
                   (is_in(n2,l) && l(n2)!=l(i,j)))
                    ImgDestination(i,j) = Imagine::WHITE;
            }
}

Imagine::Image<Imagine::Color> LoadImage(const char* img) {
    // Loading the image
    // Test to ensure the image has been loaded
    Imagine::Image<Imagine::Color> Img;
    if(!Imagine::load(Img, img)) {
        std::cout << "Image loading error!" << std::endl;
        Imagine::anyClick();
    }

    return Img;
}

void SaveImage(const Imagine::Image<Imagine::Color>& img, const char* name) {
    if(! Imagine::save(img, name)) {
        std::cerr << "Failed saving image " << name << std::endl;
        throw "Error saving image";
    }
}

/// Putting Image Img in Window W, subwindow subwin
void DisplayImage(const Imagine::Image<Imagine::Color>& Img,
                  Imagine::Window W, int subwin) {
    Imagine::setActiveWindow(W, subwin);
    display(Img);
}

void ImageSLICingAlgorithm(const Imagine::Image<Imagine::Color>& Img,
                           Imagine::Image<Imagine::Color>& ImgDestination,
                           int m, int K,
                           bool displayBorders, bool displaySuperpixels) {
    const int w=Img.width(), h=Img.height();
    Imagine::Image<int> l(w,h);

    std::vector<Superpixel> sp = SLIC(Img, l, m, K);
    enforceConnectivity(sp, l, Img);

    MakeSLICImage(displaySuperpixels, displayBorders, ImgDestination, sp, l);

    // The following calls are unused here. They are just included to
    // demonstrate their usage.
    Imagine::Coords<2>* pixels = fillSuperpixels(sp, l);
    delete [] pixels;

    std::vector< std::set<int> > adjMatrix;
    adjacencySuperpixels(l, adjMatrix);
}

int main(int argc, char* argv[]) {

    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " imgIn imgOut" << std::endl;
        return 1;
    }

    // Loading the image
    // This function raises an error if the loading fails
    Imagine::Image<Imagine::Color> Img = LoadImage(argv[1]);
    const int w=Img.width(), h=Img.height();

    // m will be the compactness parameter
    // K will be the number of superpixels
    int m, K;
    bool displayBorders;
    bool displaySuperpixels;

    GetSLICInputs(m, K, displayBorders, displaySuperpixels);

    // Preparing the destination Image
    Imagine::Image<Imagine::Color> DestinationImg = Img.clone();

    std::string names[] = {"Original image", "SLICced image"};
    Imagine::Window W = Imagine::openComplexWindow(w + 20, h + 20, "SLIC 101",
                                                   2, names);

    DisplayImage(Img, W, 0);

    ImageSLICingAlgorithm(Img, DestinationImg, m, K, displayBorders, displaySuperpixels);

    DisplayImage(DestinationImg, W, 1);

    SaveImage(DestinationImg, argv[2]);

    Imagine::endGraphics();
    return 0;
}
