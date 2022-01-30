# ClosestPointMethod
Simple C++ Implementation of the Closest Point Method for Surface PDEs

# Building
Clone using the following command to ensure polyscope is also included for visualization: 
git clone --recurse-submodules https://github.com/nathandking/ClosestPointMethod.git

Cmake is used to build the project. You should simply need to run the following commands: 
mkdir build 
cd build/ 
cmake .. 
make 

# Running examples
Once you build the examples given in the examples/ directory will be built in build/bin/. To run an example simply run the name of the example as a command. For example: 
cd build/bin/ 
PoissonCircle 
