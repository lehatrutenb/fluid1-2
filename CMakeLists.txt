cmake_minimum_required (VERSION 3.8)
project (FluidProject)

set (CMAKE_CXX_STANDARD 23)
set (CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_definitions(-DCMAKE_BUILD)

set(TYPES
    ""
    CACHE STRING "Types to precompile")
add_compile_definitions(-DTYPES=${TYPES})

set(SIZES
    ""
    CACHE STRING "Sizes to precompile")
add_compile_definitions(-DSIZES=${SIZES})

set(SAVE_DIR
    "saves/"
    CACHE STRING "dir to save ticks")
add_compile_definitions(-DSAVE_DIR="${SAVE_DIR}")

option(ToPrint "Set want to print field during work" OFF)
if (ToPrint)
    add_compile_definitions(-DTO_PRINT)
endif()

add_executable(FluidProject fluid/fluid.cpp)

option(SetDebugBuildType "Set on if want not to build as release and use O2" OFF)

if (NOT SetDebugBuildType)
    set(CMAKE_BUILD_TYPE Release)
    add_compile_options(-O2)
else()
    add_compile_definitions(fsanitize=address,undefined)
endif()

add_subdirectory(fluid/)