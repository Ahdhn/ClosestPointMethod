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

        void SetupMatrices(Grid &g, Surface &s, SpMat &E, SpMat &L);
        void GetPlotInterpMatrix(Grid &g, Surface &s, SpMat &plotE, int &Np);

        VectorXd SetInitialCondition(Surface &s);
        VectorXd SetInitialCondition3D(Surface &s);

        SpMat LaplacianSharp(SpMat &L, SpMat &E);
        SpMat ImplicitEulerMatrix(SpMat &L, SpMat &E, double &dt);
        SpMat ImplicitEulerMatrix(SpMat &L, SpMat &E, double &dt, vector<bool> &identityRows);

        double ConvergenceStudy(int nStart, int numLevels, function<vector<double>(int&)> f);

    private:
        void GetClosestPointsAndBand(Surface &s, Grid &g);

        void GetLaplacianMatrix(Grid &g, SpMat &L);
        void GetExtensionMatrix(Grid &g, Surface &s, SpMat &E);
        
        void SetIdentityRows(SpMat &A, vector<bool> &identityRows);

        vector<double> ConvergenceOrder(vector<vector<double>> &error);
        double AverageConvergenceOrder(vector<vector<double>> &error);
};
