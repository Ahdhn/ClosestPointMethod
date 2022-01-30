add_executable(HeatArc HeatArc.cpp)

target_link_libraries(HeatArc
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS HeatArc RUNTIME DESTINATION ${REL})

set_target_properties(HeatArc PROPERTIES FOLDER ${EXAMPLES_FOLDER})