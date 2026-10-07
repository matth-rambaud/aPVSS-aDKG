cmake_minimum_required (VERSION 3.12)

project (ADAPTATIVE_PVSS VERSION 0.0.1 LANGUAGES CXX)

################################################################################
# Add cmake.utils to the list of CMAKE_MODULE_PATH
################################################################################
set (MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake.utils")
list (APPEND CMAKE_MODULE_PATH "${MODULE_DIR}")

# Set C++ standard to c++17
##############################################################################
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

##############################################################################
# Print information (to ease debugging)
##############################################################################
message (STATUS "Running on system ${CMAKE_HOST_SYSTEM} with processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
if (CMAKE_CROSSCOMPILING)
  message (STATUS "Targeting system ${CMAKE_SYSTEM} with processor ${CMAKE_SYSTEM_PROCESSOR}")
endif ()
message (STATUS "Using CMake ${CMAKE_VERSION}")
message (STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}
${CMAKE_CXX_COMPILER_VERSION}")

set (CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -DBICYCL_WITH_PTHREADS")
  
################################################################################
# Looking for GMP, OpenSSL, BICYCL, Relic
################################################################################
find_package (gmp REQUIRED)
find_package (OpenSSL REQUIRED)
find_package (BICYCL REQUIRED)

add_subdirectory (benchs)

