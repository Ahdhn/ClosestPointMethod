#pragma once

#include <vector>
#include <iostream>
#include <fstream>

#include "Grid.h"
#include "Surface.h"
#include "Interpolation.h"
#include "Helpers.h"

using namespace std;
using namespace Eigen;

class Tests {
    public:
        Tests();
    
        Helpers h;

        void StandardInterpolationOrder();
        vector<double> InterpolationOrder(int n, int interpPolyDegree);
        vector<double> InterpolationOrder3D(int n, int interpPolyDegree);
};
