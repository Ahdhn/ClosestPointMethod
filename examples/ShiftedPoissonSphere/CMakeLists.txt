add_executable(ShiftedPoissonSphere ShiftedPoissonSphere.cpp)

target_link_libraries(ShiftedPoissonSphere
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ShiftedPoissonSphere RUNTIME DESTINATION ${REL})

set_target_properties(ShiftedPoissonSphere PROPERTIES FOLDER ${EXAMPLES_FOLDER})