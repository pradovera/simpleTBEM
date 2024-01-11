include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(parametrizations)

add_subdirectory(quadrature)

add_subdirectory(utilities)

add_library(single_layer STATIC single_layer.cpp)
target_link_libraries(single_layer parametrizations)
add_library(double_layer STATIC double_layer.cpp)
target_link_libraries(double_layer parametrizations)
add_library(hypersingular STATIC hypersingular.cpp)
target_link_libraries(hypersingular parametrizations)
add_library(mass_matrix STATIC mass_matrix.cpp)
target_link_libraries(mass_matrix parametrizations)
