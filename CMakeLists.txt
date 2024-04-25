cmake_minimum_required (VERSION 2.6)

# Maps to Visual Studio solution file (Tutorial.sln)
# The solution will have all targets (exe, lib, dll) 
# as Visual Studio projects (.vcproj)
PROJECT(BADChIMP)

SET(CXX_STD 17) 

#SET(LBMAIN laplace_fluid)
#SET(LBMAIN laplace_pressure)
#SET(LBMAIN six_dof)
#SET(LBMAIN co2)
#SET(LBMAIN std_case)

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )
# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#MPI STUFF
IF(WIN32)
if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
        set(MPI_MSMPI_LIB_PATH_ENV_NATIVE "$ENV{MSMPI_LIB64}")
        file(TO_CMAKE_PATH "${MPI_MSMPI_LIB_PATH_ENV_NATIVE}" MPI_MSMPI_LIB_PATH)

        # Microsoft MPI might use backslashes in the environment variables,
        # so it's important to convert to CMake-standard forward slashes
        # before appending a subdirectory with a forward slash.
        set(MPI_MSMPI_INC_PATH_ENV_NATIVE "$ENV{MSMPI_INC}")
        file(TO_CMAKE_PATH "${MPI_MSMPI_INC_PATH_ENV_NATIVE}" MPI_MSMPI_INC_PATH_ENV)
        set(MPI_MSMPI_INC_PATH_EXTRA "${MPI_MSMPI_INC_PATH_ENV}/x64")
      else()
        set(MPI_MSMPI_LIB_PATH "$ENV{MSMPI_LIB32}")

        set(MPI_MSMPI_INC_PATH_ENV_NATIVE "$ENV{MSMPI_INC}")
        file(TO_CMAKE_PATH "${MPI_MSMPI_INC_PATH_ENV_NATIVE}" MPI_MSMPI_INC_PATH_ENV)
        set(MPI_MSMPI_INC_PATH_EXTRA "${MPI_MSMPI_INC_PATH_ENV}/x86")
      endif()
ENDIF(WIN32)

find_package(MPI)
IF(MPI_CXX_FOUND)
        INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
        LIST(APPEND SCR_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES})
ENDIF(MPI_CXX_FOUND)



#if ( CMAKE_COMPILER_IS_GNUCC )
    # This is (not) tested
    # set_property( TARGET LBSOLVER APPEND_STRING PROPERTY COMPILE_FLAGS
    # -Wwrite-strings )
    #  add_definitions("-pg")
#endif ( CMAKE_COMPILER_IS_GNUCC )
#######################################################
## For debugging and memory checking on Linux        ##
## run with >valgrind --leak-check=full <executable> ##
#######################################################
IF(WIN32)
# Set compiler flags and options. 
# Here it is setting the Visual Studio warning level to 4
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 ")
ELSEIF(UNIX AND NOT APPLE)
  IF(CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -pg -std=c++${CXX_STD} -Wall")
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -std=c++${CXX_STD} -Wall")
  ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -pg -std=c++11")
  #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -std=c++${CXX_STD}")
  #SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -O3 -std=c++${CXX_STD}")
ELSE()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++${CXX_STD} -Wall")
ENDIF(WIN32)

# Command to output information to the console
# Useful for displaying errors, warnings, and debugging

message ("cxx Flags: " ${CMAKE_CXX_FLAGS})

# Sub-directories where more CMakeLists.txt exist
ADD_SUBDIRECTORY( src/lbsolver/ )
ADD_SUBDIRECTORY( src/io/       )


if(EXISTS "${PROJECT_SOURCE_DIR}/src/${LBMAIN}")
 string(COMPARE EQUAL "${LBMAIN}" "" result)
 if(result)
  message("CMAKE: NO ARGUMENTS GIVEN")
  message("CMAKE: ADDING MAIN FILE STD_CASE")
  ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/src/std_case/  )
 else()
  message("CMAKE: ADDING MAIN FILE ${LBMAIN}")
  ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/${LBMAIN}  )
 endif()
else()
 message("CMAKE: COULD NOT FIND FOLDER: ${LBMAIN}")
 message("CMAKE: ADDING MAIN FILE STD_CASE")
 ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/src/std_case/  )
endif()


message("BUILD:" ${CMAKE_BUILD_TYPE})





