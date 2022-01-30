#include <vector>
#include <iostream>

#include "Grid.h"
#include "Surface.h"
#include "Helpers.h"

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"

using namespace std;
using namespace Eigen;

Helpers h;

vector<double> ShiftedPoissonCircle(int n, bool visualize = true)
{
    Grid g(3 /* interpolation order */, 2 /* embedding space dimension */);
    Surface s("Circle");
    SpMat L, E;
    
    h.SetupMatrices(n, g, s, E, L);

    // function in the embedding space
    vector<double> ug;
    double theta;
    for(int i = 0; i < g.sizeBand; ++i)
    {
        theta = atan2(s.cpx[i][1], s.cpx[i][0]);
        ug.push_back(2 * sin(theta) + 145 * sin(12 * theta));
    }

    Map<VectorXd> b(ug.data(), ug.size());
    ug.clear();

    ///////////////////////////////////////////
    // solve the Poisson equation on a circle
    ///////////////////////////////////////////
    SpMat M = h.LaplacianSharp(L, E);
    SpMat I(M.rows(), M.cols());
    I.setIdentity();
    SpMat A = I-M;

    SparseLU<SpMat> solver;
    solver.compute(A);
    if(solver.info()!=Success) {
        cout << "decomposition failed" << endl;
    }

    VectorXd u = solver.solve(b);
    if(solver.info()!=Success) {
        cout << "solver failed" << endl;
    }

    // compute error with exact solution on the surface
    int Np = 1000;
    SpMat plotE(Np, g.sizeBand);
    h.GetPlotInterpMatrix(g, s, plotE, Np);
    VectorXd uplot = plotE * u;
    VectorXd uexact(Np);
    for(int i = 0; i < Np; ++i)
    {
        uexact[i] = sin(s.thetap[i]) + sin(12 * s.thetap[i]);
    }
    VectorXd error = uplot - uexact;
    
    // Visualize with polyscope
    if(visualize)
    {
        string curveName = "Surface dx=" + to_string(g.dx[0]);
        polyscope::registerCurveNetworkLoop2D(curveName, s.xp);
        polyscope::getCurveNetwork(curveName)->addNodeScalarQuantity("uplot", uplot);
        polyscope::getCurveNetwork(curveName)->addNodeScalarQuantity("uexact", uexact);
        polyscope::getCurveNetwork(curveName)->addNodeScalarQuantity("error", error);
    }
    
    vector<double> error_info(2);
    error_info[0] = error.lpNorm<Infinity>();
    error_info[1] = g.dx[0];
    return error_info;
}

int main(int argc, char** argv)
{
    bool visualize = true; // visualize the solution using polyscope
    if(visualize){ polyscope::init(); }

    // ShiftedPoissonCircle(6, visualize);  // uncomment if you want to run just a single example with a single dx value 

    // make a lambda for convergence study, we need function format double f(double)
    auto Poisson = [&visualize](int& n)
    {
        return ShiftedPoissonCircle(n, visualize);
    };

    h.ConvergenceStudy(6, 5, Poisson);

    if(visualize){ polyscope::show(); }
}
