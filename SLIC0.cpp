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


void ImageSLICingAlgorithm(Imagine::Image<Imagine::Color> Img, Imagine::Image<Imagine::Color> ImgDestination, int m, int K, int w, int h, bool displayBorders, bool displaySuperpixels) {

    Imagine::Image<int> l(w, h);

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                  SLIC algorithm                                 ////
    //// ------------------------------------- </> ------------------------------------- ////

    std::vector<Superpixel> Superpixels = SLIC(Img, l, m, K, h, w);

    //// ------------------------------------- </> ------------------------------------- ////
    ////                              Ensuring connectivity                              ////
    //// ------------------------------------- </> ------------------------------------- ////

    ConnectivityEnforcement(K, w, h, l, Superpixels, Img);

    //// ------------------------------------- </> ------------------------------------- ////
    ////                               Output conditioning ?                             ////
    //// ------------------------------------- </> ------------------------------------- ////

    // Haven't determined yet what to return, maybe returning l is enough

    for(int i=0; i < w; i++) {
        for(int j = 0; j < h; j++) {
            // This loop fills each Superpixel's _contents vector with the pixels it contains, as indicated by the values in l
            //
            // Recall : l(i, j) is an int corresponding to the index of the Superpixel containing Img(i, j) in Superpixels
            // e. g. : l(i, j) = k means that Img(i, j) belongs to Superpixels[k]

            Superpixels[l(i, j)].push_contents(Imagine::Coords<2>(i, j));
        }
    }

    //// ------------------------------------- </> ------------------------------------- ////
    ////                              Making the SLIC image                              ////
    //// ------------------------------------- </> ------------------------------------- ////

    MakeSLICImage(displaySuperpixels, displayBorders, ImgDestination, Img, w, h, l, Superpixels);
    // Maybe it would be better to return l, or even the SLICed image (even if it is saved) ?
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
    Imagine::Window W = Imagine::openComplexWindow(w + 20, h + 20, "SLIC 101", 2, names);   // A window with 2 subwindows

    DisplayImage(Img, W, 0);

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                      SLIC                                       ////
    //// ------------------------------------- </> ------------------------------------- ////

    ImageSLICingAlgorithm(Img, DestinationImg, m, K, w, h, displayBorders, displaySuperpixels);

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                 Display and save                                ////
    //// ------------------------------------- </> ------------------------------------- ////

    DisplayImage(DestinationImg, W, 1);

    SaveImage(DestinationImg, argv[2]);

    Imagine::endGraphics();
    return 0;
}
