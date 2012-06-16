cmake_minimum_required(VERSION 2.6)
project(vef)
enable_language(Fortran)

#----------------------------BEGIN COMPILER STUFF-------------------------------
# uncomment next line to force Intel's compiler
#set (CMAKE_Fortran_COMPILER "ifort")
# uncomment next line to force IBM's compiler
#set (CMAKE_Fortran_COMPILER "xlf")
# uncomment next line to force GNU's compiler
set (CMAKE_Fortran_COMPILER "gfortran")
# uncomment next line to force NAG's compiler
#set (CMAKE_Fortran_COMPILER "nagfor")

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
 set (CMAKE_Fortran_FLAGS "-O0 -g -Wall -pedantic")
 set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -pedantic")
 set (CMAKE_Fortran_LINK_FLAGS "-O0 -g -Wall")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

if (Fortran_COMPILER_NAME STREQUAL "ifort")
 set (CMAKE_Fortran_FLAGS "-O0 -g -warn all")
 set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
 set (CMAKE_Fortran_LINK_FLAGS "-O0 -g -mkl")
endif (Fortran_COMPILER_NAME STREQUAL "ifort")

if (Fortran_COMPILER_NAME STREQUAL "nagfor")
 set (CMAKE_Fortran_FLAGS "-O0 -g90 -C=all")
 set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g90 -C=all")
 set (CMAKE_Fortran_LINK_FLAGS "-O0 -g90")
endif (Fortran_COMPILER_NAME STREQUAL "nagfor")

message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#-----------------------------END COMPILER STUFF--------------------------------


#---------------------------BEGIN ADD EXTRA MODULES-----------------------------
# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable, which
# stores the directories of additional CMake modules
# (e.g. MacroOutOfSourceBuild.cmake):
set(CMAKE_MODULE_PATH ${vef_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# The macro below forces the build directory to be different from source 
# directory
include(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build. Please, see README.")
#----------------------------END ADD EXTRA MODULES------------------------------



#------------------------------BEGIN ADD LIBRARIES------------------------------
find_package(BLAS)
find_package(LAPACK)
#-------------------------------END ADD LIBRARIES-------------------------------



#-----------------------------BEGIN ADD SOURCE CODE-----------------------------
add_subdirectory(src)
#------------------------------END ADD SOURCE CODE------------------------------

