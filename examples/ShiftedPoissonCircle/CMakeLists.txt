add_executable(ShiftedPoissonCircle ShiftedPoissonCircle.cpp)

target_link_libraries(ShiftedPoissonCircle
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ShiftedPoissonCircle RUNTIME DESTINATION ${REL})

set_target_properties(ShiftedPoissonCircle PROPERTIES FOLDER ${EXAMPLES_FOLDER})