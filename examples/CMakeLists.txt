include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(assembler assembler.cpp)
add_executable(mass mass.cpp)
add_executable(solver solver.cpp)

target_link_libraries(assembler parametrizations solvers)
target_link_libraries(mass parametrizations mass_matrix)
target_link_libraries(solver parametrizations solvers)

