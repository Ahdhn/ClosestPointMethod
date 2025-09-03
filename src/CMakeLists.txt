add_library(src
            FDMatrices.cpp
            Grid.cpp
            Helpers.cpp
            Interpolation.cpp
            Node.cpp
            Surface.cpp
			FDMatrices.h
            Grid.h
            Helpers.h
            Interpolation.h
            Node.h
            Surface.h)

target_link_libraries(src 
                      Eigen3::Eigen 
                      polyscope)

target_include_directories(src PUBLIC 
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include>)