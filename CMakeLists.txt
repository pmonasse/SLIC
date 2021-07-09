cmake_minimum_required(VERSION 2.6)

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
target_link_libraries(SLIC ${JPEG_LIBRARIES} ${PNG_LIBRARIES})
