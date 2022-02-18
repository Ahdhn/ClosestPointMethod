#include "Helpers.h"


Helpers::Helpers()
{

}

// construct grid, compute closest points, build discrete Laplacian L, and build closest point extension E
void Helpers::SetupMatrices(Grid &g, Surface &s, SpMat &E, SpMat &L)
{
    GetClosestPointsAndBand(s, g);

    L.resize(g.sizeBand(), g.sizeBand());
    GetLaplacianMatrix(g, L);

    E.resize(s.cpx().size(), g.sizeBand());
    GetExtensionMatrix(g, s, E);
}


// construct plotting interpolation matrix
void Helpers::GetPlotInterpMatrix(Grid &g, Surface &s, SpMat &plotE, int &Np)
{
    s.GetSurfaceParameterization(Np);

    Interpolation plotInterp(g, s.xp());
    plotInterp.BuildInterpolationMatrix(g, plotE);
}


VectorXd Helpers::SetInitialCondition(Surface &s)
{
    // function in the embedding space
    VectorXd u(s.cpx().size());
    double theta;
    for(int i = 0; i < s.cpx().size(); ++i)
    {
        theta = atan2(s.cpx()[i][1], s.cpx()[i][0]);
        u[i] = sin(theta);
    }

    return u;
}


VectorXd Helpers::SetInitialCondition3D(Surface &s)
{
    // function in the embedding space
    VectorXd u(s.cpx().size());
    double theta;
    double phi;
    for(int i = 0; i < s.cpx().size(); ++i)
    {
        theta = atan2(s.cpx()[i][1], s.cpx()[i][0]);
        phi = acos(s.cpx()[i][2] / s.surfaceParams()[0]);
        u[i] = sin(theta)*sin(phi);
    }

    return u;
}


// Numerically stable discrete Laplace-Beltrami operator, see Macdonald & Ruuth 2009 (section 2.2.3 of https://steveruuth.org/wp-content/uploads/2020/10/icpm.pdf)
SpMat Helpers::LaplacianSharp(SpMat &L, SpMat &E)
{
    SpMat M;

    SpMat LnoDiag = L;
    LnoDiag -= L.diagonal().asDiagonal();

    SpMat LnoDiagE;
    LnoDiagE = LnoDiag * E;
    SpMat Ldiag;
    Ldiag = L.diagonal().asDiagonal();
    Ldiag += LnoDiagE;
    M = Ldiag;

    return M;
}


// Implicit Euler time-stepping matrix
SpMat Helpers::ImplicitEulerMatrix(SpMat &L, SpMat &E, double &dt)
{
    SpMat M = LaplacianSharp(L, E);
    SpMat I(M.rows(), M.cols());
    I.setIdentity();

    return I - dt * M;
}


// Implicit Euler time-stepping matrix with some identity rows
SpMat Helpers::ImplicitEulerMatrix(SpMat &L, SpMat &E, double &dt, vector<bool> &identityRows)
{
    SpMat A = ImplicitEulerMatrix(L, E, dt);
    SetIdentityRows(A, identityRows);

    return A;
}


// run a convergence study for a function f that returns an error value for all dx specified
// nStart: coarsest grid resolution power
// numLevels: number different grid resolutions
// f: any function that returns an error and a given grid resolution
double Helpers::ConvergenceStudy(int nStart, int numLevels, function<vector<double>(int&)> f)
{
    // compute all grid division powers for numLevels
    vector<int> gridDivisionPower(numLevels);
    for(int i = 0; i < numLevels; ++i)
    {
        gridDivisionPower[i] = i + nStart;
    }

    // compute error for each dx value for the function f
    vector<vector<double>> error(numLevels, vector<double>(2));
    for(int i = 0; i < numLevels; ++i)
    {
        error[i] = f(gridDivisionPower[i]); // must return error and dx 
    }
    
    vector<double> order = ConvergenceOrder(error);

    // output table of results
    cout << setw(8) << "dx" << '\t' << setw(8) << "error" << '\t' << setw(8) << "order" << endl;
    for(int i = 0; i < numLevels; ++i)
    {
        cout << setw(8) << error[i][1] << '\t' << setw(8) << error[i][0] << '\t' << setw(8) << order[i] << endl;
    }

    double avg_conv_order = AverageConvergenceOrder(error);
    cout << "Average Convergence Order = " << avg_conv_order << endl;

    return avg_conv_order;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//          Private Functions
/////////////////////////////////////////////////////////////////////////////////////////////////


void Helpers::GetClosestPointsAndBand(Surface &s, Grid &g)
{
    // compute closest points and distances to all grid points    
    s.GetClosestPoints(g.xg());

    // compute the banded domain based on the distance function
    g.SetComputationalBand(s.dist());
    s.BandClosestPoints(g);
}


// build Laplacian matrix for grid points within band
void Helpers::GetLaplacianMatrix(Grid &g, SpMat &L)
{
    FDMatrices mFDMat;
    mFDMat.BuildLaplacianMatrix(g, L);
}


// interpolation on the banded grid
void Helpers::GetExtensionMatrix(Grid &g, Surface &s, SpMat &E)
{
    Interpolation interp(g, s.cpx());
    interp.BuildInterpolationMatrix(g, E);
}


// Make rows in matrix rows of the identity matrix when identityRows[i] == true
void Helpers::SetIdentityRows(SpMat &A, vector<bool> &identityRows)
{
    for (int k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            if(identityRows[it.row()])
            {
                if(it.row() == it.col())
                {
                    it.valueRef() = 1.0;
                }
                else
                {
                    it.valueRef() = 0.0;
                }
            }
        }
    }
}


// Compute the order of convergence
vector<double> Helpers::ConvergenceOrder(vector<vector<double>> &error)
{
    vector<double> order(error.size());
    order[0] = NAN;
    for(int i = 1; i < error.size(); ++i)
    {
        order[i] = log(error[i-1][0]/error[i][0]) / log(error[i-1][1]/error[i][1]);
    }

    return order;
}


double Helpers::AverageConvergenceOrder(vector<vector<double>> &error)
{
    vector<double> order = ConvergenceOrder(error);
    return accumulate(order.begin() + 1, order.end(), 0.0) / (order.size() - 1);
}
