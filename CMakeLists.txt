cmake_minimum_required(VERSION 2.6)
file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED)

project(SLIC)

find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

add_executable(SLIC
  image.h
  io_jpg.h io_jpg.c
  io_png.h io_png.c
  slic.h slic.cpp
  main.cpp cmdLine.h)
ImagineUseModules(SLIC Images)
target_link_libraries(SLIC ${JPEG_LIBRARIES} ${PNG_LIBRARIES})
