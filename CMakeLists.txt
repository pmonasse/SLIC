cmake_minimum_required(VERSION 2.6)
file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED)

project(SLIC-v0)

add_executable(SLIC-v0 
               superpixel.h
               superpixel.cpp
               SLICTools.cpp
               SLICTools.h
               SLIC0.cpp
               )

ImagineUseModules(SLIC-v0 Images)
