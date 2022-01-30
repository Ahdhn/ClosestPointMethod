#pragma once

#include <vector>
#include <iostream>
#include <iomanip> // for setw
#include <numeric> // for accumulate

#include "Grid.h"
#include "Surface.h"
#include "FDMatrices.h"
#include "Interpolation.h"

using namespace std;
using namespace Eigen;

class Helpers {
    public:
        Helpers();

        void SetOriginalGrid(Grid &g, double (&bBox)[2], int n);
        void GetClosestPointsAndBand(Surface &s, Grid &g);
    
        VectorXd SetInitialCondition(Surface &s);
        VectorXd SetInitialCondition3D(Surface &s);

        void GetLaplacianMatrix(Grid &g, SpMat &L);
        void GetExtensionMatrix(Grid &g, Surface &s, SpMat &E);
        void GetPlotInterpMatrix(Grid &g, Surface &s, SpMat &plotE, int &Np);
    
        void SetupMatrices(int n, Grid &g, Surface &s, SpMat &E, SpMat &L);

        SpMat LaplacianSharp(SpMat &L, SpMat &E);
        SpMat ImplicitEulerMatrix(SpMat &L, SpMat &E, double &dt);
        SpMat ImplicitEulerMatrix(SpMat &L, SpMat &E, double &dt, vector<bool> &identityRows);
        void SetIdentityRows(SpMat &A, vector<bool> &identityRows);

        void ConvergenceStudy(int nStart, int numLevels, function<vector<double>(int&)> f);
        vector<double> ConvergenceOrder(vector<vector<double>> &error);
        double AverageConvergenceOrder(vector<vector<double>> &error);
};
