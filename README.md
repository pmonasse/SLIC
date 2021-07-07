Superpixels through K-means (SLIC)
by Robin Gay <robin.gay@eleves.enpc.fr>
and Pascal Monasse <pascal.monasse@enpc.fr>

Build
-----
Prerequisites: CMake (https://cmake.org/), Imagine++ (http://imagine.enpc.fr/~monasse/Imagine++/)

$ mkdir Build && cd Build && cmake -DCMAKE_BUILD_TYPE:bool=Release ..
$ cmake --build .

Usage
-----
./SLIC [options] imgIn imgOut
-k ARG Required number of superpixels (1000)
-m ARG Compactness parameter (100)
