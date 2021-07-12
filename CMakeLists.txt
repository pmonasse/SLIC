cmake_minimum_required(VERSION 2.6)

project(SLIC)

# zlib is needed by libPNG
find_package(ZLIB)
if(NOT ZLIB_FOUND)
  add_subdirectory(third_party/zlib)
  set(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
endif(NOT ZLIB_FOUND)
include_directories(${ZLIB_INCLUDE_DIRS})

find_package(PNG)
if(PNG_FOUND)
  add_definitions(${PNG_DEFINITIONS})
else(PNG_FOUND)
  set(SKIP_INSTALL_ALL 1)
  add_subdirectory(third_party/png)
endif(PNG_FOUND)
include_directories(${PNG_INCLUDE_DIRS})

find_package(JPEG)
if(NOT JPEG_FOUND)
  add_subdirectory(third_party/jpeg)
endif(NOT JPEG_FOUND)
include_directories(${JPEG_INCLUDE_DIR})

add_executable(SLIC
  image.h
  io_jpg.h io_jpg.c
  io_png.h io_png.c
  slic.h slic.cpp
  main.cpp cmdLine.h)
target_link_libraries(SLIC ${JPEG_LIBRARIES} ${PNG_LIBRARIES})
