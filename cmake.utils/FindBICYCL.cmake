# FindBICYCL
# -----------
#
# The module defines the following variables:
#
#     BICYCL_FOUND
#     BICYCL_INCLUDE_DIRS
#
# and the following imported target (if it does not already exist):
#
#  BICYCL::BICYCL - The BICYCL library
#
#
# Requires CMake >= 3.0

set (BICYCL_DIR "${BICYCL_DIR}" CACHE PATH "Directory to search for BICYCL")

# Look for the include directory
find_path (BICYCL_INC_DIR NAMES bicycl.hpp
                         HINTS "${_libdir}/.." "${BICYCL_DIR}"
                         PATH_SUFFIXES include)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (BICYCL DEFAULT_MSG BICYCL_INC_DIR)

if (BICYCL_FOUND)
  set (BICYCL_INCLUDE_DIRS "${BICYCL_INC_DIR}")
  mark_as_advanced (BICYCL_DIR)
  if (NOT TARGET BICYCL::BICYCL)
    # For now we make the target global, because this file is included from a
    # CMakeLists.txt file in a subdirectory. With CMake >= 3.11, we could make
    # it global afterwards with
    # set_target_properties(BICYCL::BICYCL PROPERTIES IMPORTED_GLOBAL TRUE)
    add_library (BICYCL::BICYCL INTERFACE IMPORTED GLOBAL)
    set_target_properties (BICYCL::BICYCL PROPERTIES
              INTERFACE_INCLUDE_DIRECTORIES "${BICYCL_INCLUDE_DIRS}" )
  endif()
endif()

mark_as_advanced(BICYCL_INC_DIR BICYCL_LIBRARY)
