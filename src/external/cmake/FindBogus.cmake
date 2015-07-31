# Find Bogus includes and libraries.
# The following variables are set if Bogus is found.  If Bogus is not
# found, Bogus_FOUND is set to false.
#  Bogus_FOUND        - True when the Bogus include directory is found.
#  Bogus_INCLUDE_DIRS - the path to where the Siconos Bogus include files are.
#  Bogus_LIBRARY_DIRS - The path to where the Siconos library files are.
#  Bogus_LIBRARIES    - The libraries to link against Siconos Bogus

# One may want to use a specific Bogus Library by setting
# Bogus_LIBRARY_DIRECTORY before FIND_PACKAGE(Bogus)
INCLUDE(FindPackageHandleStandardArgs)

IF(Bogus_LIBRARY_DIRECTORY)
  FIND_LIBRARY(Bogus_LIBRARY bogus PATHS "${Bogus_LIBRARY_DIRECTORY}")
ELSE(Bogus_LIBRARY_DIRECTORY)
  FIND_LIBRARY(Bogus_LIBRARY bogus)
ENDIF(Bogus_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Bogus
  REQUIRED_VARS Bogus_LIBRARY)

IF(Bogus_LIBRARY)
  GET_FILENAME_COMPONENT(Bogus_LIBRARY_DIRS ${Bogus_LIBRARY} PATH)
  SET(Bogus_LIBRARIES ${Bogus_LIBRARY} ${Bogus_COMMON_LIBRARY})
  GET_FILENAME_COMPONENT(Bogus_LIBRARY_DIRS_DIR ${Bogus_LIBRARY_DIRS} PATH)

  FIND_PATH(Bogus_INCLUDE_DIRS Core/Block.hpp
    HINTS ${Bogus_LIBRARY_DIRS_DIR} ${Bogus_LIBRARY_DIRS_DIR_DIR} 
    ENV PATH
    PATH_SUFFIXES include/bogus)

ELSE(Bogus_LIBRARY)
  IF(Bogus_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Bogus library not found. Please specify library location in Bogus_LIBRARY_DIRECTORY")
  ENDIF(Bogus_FIND_REQUIRED)
ENDIF(Bogus_LIBRARY)
