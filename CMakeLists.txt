cmake_minimum_required(VERSION 3.10)

# Set the project name
project(LinearLeastSquares)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(LinearLeastSquares main.cpp matrix.cpp linear_least_squares.cpp eigen.cpp nonlinear.cpp)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Optionally link against libraries (if you need any)
# target_link_libraries(LinearLeastSquares <library_name>)

