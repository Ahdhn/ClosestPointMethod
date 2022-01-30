#pragma once

#include <Eigen/Sparse>
#include "Grid.h"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

using namespace std;

class FDMatrices{
    public:
        FDMatrices();

        void BuildLaplacianMatrix(Grid &g, SpMat &L);
        void BuildGradientMatrices(Grid &g, vector<SpMat> &Dc);
};
