cmake_minimum_required(VERSION 2.6)

project(SLIC)

find_package(PNG)
if(NOT PNG_FOUND)
  # zlib is needed by libPNG
  find_package(ZLIB)
  if(NOT ZLIB_FOUND)
    add_subdirectory(third_party/zlib-1.2.13)
    set(ZLIB_FOUND TRUE)
  endif()
  include_directories(${ZLIB_INCLUDE_DIRS})
  set(SKIP_INSTALL_ALL 1)
  add_subdirectory(third_party/libpng-1.6.39)
endif()
add_definitions(${PNG_DEFINITIONS})
include_directories(${PNG_INCLUDE_DIRS})

find_package(JPEG)
if(NOT JPEG_FOUND)
  add_subdirectory(third_party/jpeg-9e)
endif()
include_directories(${JPEG_INCLUDE_DIR})

add_executable(SLIC
  image.h
  io_jpg.h io_jpg.c
  io_png.h io_png.c
  slic.h slic.cpp
  main.cpp cmdLine.h)
target_link_libraries(SLIC ${JPEG_LIBRARIES} ${PNG_LIBRARIES})
