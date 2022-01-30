add_executable(HeatSphere HeatSphere.cpp)

target_link_libraries(HeatSphere
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS HeatSphere RUNTIME DESTINATION ${REL})

set_target_properties(HeatSphere PROPERTIES FOLDER ${EXAMPLES_FOLDER})