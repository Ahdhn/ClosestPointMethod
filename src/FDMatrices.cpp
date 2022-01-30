#include "FDMatrices.h"

FDMatrices::FDMatrices()
{

}


// TO DO: make this work for different grid spacing in each dimension, but would that be less efficient?
void FDMatrices::BuildLaplacianMatrix(Grid &g, SpMat &L)
{
    vector<T> coeffs; // list of non-zeros coefficients, (row, col, value) triplets to be used to construct sparse matrix afterward
    for(int row = 0; row < g.sizeBand; ++row)
    {
        coeffs.push_back(T(row, row, -2.0 * (double)g.dim)); // diagonal entry
        for(int nbr = 0; nbr < 2 * g.dim; ++nbr)
        {
            if(g.bandNodes[row].neighbours[nbr] != nullptr)
            {
                coeffs.push_back(T(row, g.bandNodes[row].neighbours[nbr]->bandIndex, 1.0)); // neighour entries
            }
        }
    }
    
    L.setFromTriplets(coeffs.begin(), coeffs.end());
    L /= (g.dx[0] * g.dx[0]);
}


void FDMatrices::BuildGradientMatrices(Grid &g, vector<SpMat> &Dc)
{
    vector<vector<T>> coeffs(g.dim); // list of non-zeros coefficients for each dimension in a vector, (row, col, value) triplets to be used to construct sparse matrix afterward

    for(int row = 0; row < g.sizeBand; ++row)
    {
        for(int nbr = 0; nbr < 2 * g.dim; ++nbr)
        {
            if(g.bandNodes[row].neighbours[nbr] != nullptr)
            {
                coeffs[nbr/2].push_back(T(row, g.bandNodes[row].neighbours[nbr]->bandIndex, pow(-1, nbr) * 0.5)); // neighour entry
            }
        }
    }
    
    for(int d = 0; d < g.dim; ++d)
    {
        Dc[d].setFromTriplets(coeffs[d].begin(), coeffs[d].end());
        Dc[d] /= g.dx[d];
    }
}
