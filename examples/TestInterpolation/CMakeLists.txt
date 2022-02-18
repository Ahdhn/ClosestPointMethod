add_executable(TestInterpolation TestInterpolation.cpp)

target_link_libraries(TestInterpolation
                      src
                      polyscope
                      Eigen3::Eigen)
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS TestInterpolation RUNTIME DESTINATION ${REL})

set_target_properties(TestInterpolation PROPERTIES FOLDER ${EXAMPLES_FOLDER})