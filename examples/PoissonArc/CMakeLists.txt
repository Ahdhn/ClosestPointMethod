add_executable(PoissonArc PoissonArc.cpp)

target_link_libraries(PoissonArc
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS PoissonArc RUNTIME DESTINATION ${REL})

set_target_properties(PoissonArc PROPERTIES FOLDER ${EXAMPLES_FOLDER})