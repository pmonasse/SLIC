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

    // Chronometer for test version
    auto time_before_CE = std::chrono::high_resolution_clock::now();
    ConnectivityEnforcement(K, w, h, l, Superpixels, Img);
    // Chronometer for test version
    auto time_after_CE = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time_CE = (time_after_CE - time_before_CE);
    std::cout << "Elapsed time for connectivity enforcement: " << elapsed_time_CE.count() << std::endl;

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

int main() {

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                     Inputs                                      ////
    //// ------------------------------------- </> ------------------------------------- ////

    // This rand initializer is useful to pick a random image from the Default_Images directory if no image is specified
    RandInit();

    std::string path;
    std::cout << "Path to the image to SLIC: ";
    std::getline(std::cin, path);

    if(path.empty()) {
        path.append(srcPath("./Default_Images"));
    }

    std::string imageName;
    std::cout << "Image name: ";
    std::getline(std::cin, imageName);

    if(imageName.empty()) {
        imageName.append(std::to_string(int(rand()) % 100));
        imageName.append(".png");
        std::cout << "Random image seleced: " << imageName << std::endl;
    }

    std::string Image;

    if(path[path.size() - 1] == '/') {
        Image = path + imageName;
    }
    else {
        path.append("/");
        Image = path + imageName;
    }

    const char* image = Image.c_str();

    // These int will be the width and height of the window
    int w, h;

    // Loading the image
    // This function raises an error if the loading fails
    Imagine::Image<Imagine::Color> Img = LoadImage(image, w, h);

    // m will be the compactness parameter
    // K will be the number of superpixels
    int m, K;
    bool displayBorders;
    bool displaySuperpixels;

    GetSLICInputs(m, K, displayBorders, displaySuperpixels);

    // Preparing the destination Image
    Imagine::Image<Imagine::Color> DestinationImg = Img.clone();

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                     Display                                     ////
    //// ------------------------------------- </> ------------------------------------- ////

    std::string names[] = {"Original image", "SLICced image"};
    Imagine::Window W = Imagine::openComplexWindow(w + 20, h + 20, "SLIC 101", 2, names);   // A window with 2 subwindows

    DisplayImage(Img, W, 0, w, h);

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                      SLIC                                       ////
    //// ------------------------------------- </> ------------------------------------- ////

    auto time_before = std::chrono::high_resolution_clock::now();
    ImageSLICingAlgorithm(Img, DestinationImg, m, K, w, h, displayBorders, displaySuperpixels);
    auto time_after = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time = (time_after - time_before);

    std::cout << "Elapsed time for SLICing: " << elapsed_time.count() << " for a " << Img.width()*Img.height() << "-pixel image" << std::endl;

    //// ------------------------------------- </> ------------------------------------- ////
    ////                                 Display and save                                ////
    //// ------------------------------------- </> ------------------------------------- ////

    DisplayImage(DestinationImg, W, 1, w, h);

    SaveSLICImage(DestinationImg, w, h, imageName);

    delete[] image;
    delete[] Img.data();
    delete[] DestinationImg.data();

    Imagine::endGraphics();
    return 0;
}
