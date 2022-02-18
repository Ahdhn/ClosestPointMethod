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

vector<double> PoissonArc(int n, bool visualize = true)
{
    Surface s("Arc");
    Grid g(3 /* interpolation order */, 2 /* embedding space dimension */, s.bBox, n);

    SpMat L, E;    
    h.SetupMatrices(g, s, E, L);

    VectorXd b = h.SetInitialCondition(s);
    
    SpMat A = h.LaplacianSharp(L, E);
    A = -A;

    // make rows corresponding to the boundary grid nodes as identity rows
    for (int k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            if(s.bdy()[it.row()] > 0)
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
    int Np = 100;
    SpMat plotE(Np, g.sizeBand());
    h.GetPlotInterpMatrix(g, s, plotE, Np);
    VectorXd uplot = plotE * u;
    vector<double> thetap_cpy = s.thetap();
    Map<VectorXd> uexact(thetap_cpy.data(), thetap_cpy.size());
    uexact = uexact.array().sin();
    VectorXd error = uplot - uexact;
    
    // Visualize with polyscope
    if(visualize)
    {
        string curveName = "Surface dx=" + to_string(g.dx()[0]);
        polyscope::registerCurveNetworkLine2D(curveName, s.xp());
        polyscope::getCurveNetwork(curveName)->addNodeScalarQuantity("uplot", uplot);
        polyscope::getCurveNetwork(curveName)->addNodeScalarQuantity("uexact", uexact);
        polyscope::getCurveNetwork(curveName)->addNodeScalarQuantity("error", error);
    }

    vector<double> error_info(2);
    error_info[0] = error.lpNorm<Infinity>();
    error_info[1] = g.dx()[0];
    return error_info;
}


int main(int argc, char** argv)
{
    bool visualize = true; // visualize the solution using polyscope
    if(visualize){ polyscope::init(); }

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    // PoissonEquationArc(6, visualize);  // uncomment if you want to run just a single example with a single dx value 

    // make a lambda for convergence study, we need function format double f(double)
    auto Poisson = [&visualize](int& n)
    {
        return PoissonArc(n, visualize);
    };

    double avg_conv_order = h.ConvergenceStudy(6, 5, Poisson);
    if(abs(avg_conv_order - 0.8016166844456067) < 1e-16)
    {
        cout << "PASSED: average convergence order test" << endl;
    }
    else
    {
        cout << "FAILED: average convergence order test" << endl;
    }

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    if(visualize){ polyscope::show(); }
}
