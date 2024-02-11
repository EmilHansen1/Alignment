cmake_minimum_required(VERSION 3.0.0)
project(rot_tdse_solver VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(rot_tdse_solver rot_tdse_solver.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
