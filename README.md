Bilateral K-Means for Superpixel Computation (the SLIC Method)
version 1.0 (2022)
by Robin Gay <robin.gay@eleves.enpc.fr>
and Pascal Monasse <pascal.monasse@enpc.fr>

Updates and future versions
---------------------------
https://github.com/pmonasse/SLIC

Build
-----
Prerequisites: CMake (https://cmake.org/)
```
$ mkdir Build && cd Build && cmake -DCMAKE_BUILD_TYPE:string=Release ..
$ cmake --build .
```

Usage
-----
```
./SLIC [options] imgIn imgOut
-k ARG Required number of superpixels (1000)
-m ARG Compactness parameter (40)
-g ARG Radius for minimal gradient search (0)
-c ARG Color of boundary (255,255,255)
```

Example
-------
```
./SLIC data/exampleIn.jpg out.jpg
'''
Compare resulting image out.jpg with data/exampleOut.jpg.

Files
-----
AUTHORS.txt     cmdLine.h  io_jpg.h     LICENSE.txt  slic.cpp(*)   third_party
image.h         io_png.c   main.cpp(*)  slic.h(*)
CMakeLists.txt  io_jpg.c   io_png.h     README.md
(*) reviewed in IPOL.
