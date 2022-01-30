#pragma once

#include <iostream>
#include <vector>
#include "math.h"
#include "Grid.h"
#include "Surface.h"

#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

using namespace std;

class Interpolation{
    public:
        Interpolation(Grid &g, vector<vector<double>> &xquery);

        int Nq;
        vector<vector<double>> xq;
        vector<int> queryWhichGrid;
        vector<vector< vector<double>>> w; // interpolation weights for each dimension
        vector<int> Ibpt;

        // computing weights
        void InterpolationWeights1D(vector<double> &x, vector<double> &w);
        void FindInterpolationBasePoint(Grid &g);

        void BuildInterpolationWeights1D(Grid &g, vector<double> &x, double &xq_subset, vector<double> &w);
        void BuildInterpolationWeights(Grid &g, vector<vector<vector<double>>> &wx);
        void BuildInterpolationMatrix(Grid &g, SpMat &E);

        // testing
        void TestRowSumOne(SpMat &E);
        void VisualizeInterpolationStencil(int &q, Grid &g, int &Ibpt, vector<double> &xq, int &queryWG);
};
