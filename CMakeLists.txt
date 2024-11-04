cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C CXX)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 17)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

message(STATUS "Building ${SKBUILD_PROJECT_NAME} for ${CMAKE_SYSTEM_PROCESSOR}")

# --- Detect C++ flags ---------------------------------------------------------

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-funroll-loops HAVE_UNROLL_LOOPS)
check_cxx_compiler_flag(-faligned-new HAVE_ALIGNED_NEW)
check_cxx_compiler_flag(-Wno-narrowing HAVE_NARROWING)
if(HAVE_UNROLL_LOOPS)
    add_compile_options(-funroll-loops)
endif()
if(HAVE_ALIGNED_NEW)
    add_compile_options(-faligned-new)
endif()
if(HAVE_NARROWING)
    add_compile_options(-Wno-narrowing)
endif()


# --- Setup include folders ----------------------------------------------------

set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Add scripts for Cython ---------------------------------------------------

include("src/scripts/cmake/CythonExtension.cmake")

# --- Compile code -------------------------------------------------------------

add_subdirectory("src")
