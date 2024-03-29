CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (vef)
ENABLE_LANGUAGE (Fortran)

#----------------------------BEGIN COMPILER STUFF-------------------------------
set (CMAKE_BUILD_TYPE "Debug")
set (CMAKE_VERBOSE_MAKEFILE false)

GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

MESSAGE ("CMAKE_Fortran_COMPILER full path: ${CMAKE_Fortran_COMPILER}")
MESSAGE ("Fortran compiler: ${Fortran_COMPILER_NAME}")
#-----------------------------END COMPILER STUFF--------------------------------


#---------------------------BEGIN ADD EXTRA MODULES-----------------------------
# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable, which
# stores the directories of additional CMake modules
# (e.g. MacroOutOfSourceBuild.cmake):
SET (CMAKE_MODULE_PATH ${vef_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# The macro below forces the build directory to be different from source
# directory
INCLUDE (MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build. Please, see README.")
#----------------------------END ADD EXTRA MODULES------------------------------


#-------------------------------BEGIN ADD OPENMP--------------------------------
FIND_PACKAGE (OpenMP)
IF (OPENMP_FOUND)
  SET (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_C_FLAGS}")
ENDIF (OPENMP_FOUND)
#--------------------------------END ADD OPENMP---------------------------------



#-----------------------------BEGIN ADD SOURCE CODE-----------------------------
ADD_SUBDIRECTORY (src)
#------------------------------END ADD SOURCE CODE------------------------------
