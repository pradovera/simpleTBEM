include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(solvers STATIC solvers.cpp)
target_link_libraries(solvers quadrature single_layer double_layer hypersingular mass_matrix complex_bessel)

