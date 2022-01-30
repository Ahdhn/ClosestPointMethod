add_executable(HeatCircle HeatCircle.cpp)

target_link_libraries(HeatCircle
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS HeatCircle RUNTIME DESTINATION ${REL})

set_target_properties(HeatCircle PROPERTIES FOLDER ${EXAMPLES_FOLDER})