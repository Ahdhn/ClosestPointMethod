#include <vector>
#include <iostream>

#include "Grid.h"
#include "Surface.h"
#include "Helpers.h"

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"

using namespace std;
using namespace Eigen;

Helpers h;

vector<double> ShiftedPoissonSphere(int n, bool visualize = true)
{
    Surface s("Sphere");
    Grid g(3 /* interpolation order */, 3 /* embedding space dimension */, s.bBox, n);
    SpMat L, E;

    h.SetupMatrices(g, s, E, L);

    // function in the embedding space
    vector<double> ug;
    double theta;
    double phi;
    for(int i = 0; i < g.sizeBand(); ++i)
    {
        theta = atan2(s.cpx()[i][1], s.cpx()[i][0]);
        if(theta < 0)
        {
            theta += 2 * M_PI;
        }
        phi = acos(s.cpx()[i][2] / s.surfaceParams()[0]);
        ug.push_back( -29 * (cos(3*theta) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1) ));
    }

    Map<VectorXd> b(ug.data(), ug.size());
    ug.clear();

    ///////////////////////////////////////////
    // solve the Poisson equation on a circle
    ///////////////////////////////////////////
    SpMat M = h.LaplacianSharp(L, E);
    SpMat I(M.rows(), M.cols());
    I.setIdentity();
    SpMat A = I+M;

    // cout << A.rows() << endl;
    // cout << A.cols() << endl;

    BiCGSTAB<SpMat> solver;
    solver.compute(A);
    if(solver.info()!=Success) {
        cout << "decomposition failed" << endl;
    }

    VectorXd u = solver.solve(b);
    if(solver.info()!=Success) {
        cout << "solver failed" << endl;
    }

    // compute error with exact solution on the surface
    int Np = 10000; // NOTE: must be some number squared
    SpMat plotE(Np, g.sizeBand());
    h.GetPlotInterpMatrix(g, s, plotE, Np);
    VectorXd uplot = plotE * u;
    VectorXd uexact(Np);
    for(int i = 0; i < Np; ++i)
    {
        uexact[i] = cos(3*s.thetap()[i]) * pow(sin(s.phip()[i]), 3) * (9 * pow(cos(s.phip()[i]),2) - 1) ;
    }
    VectorXd error = uplot - uexact;

    // Visualize with polyscope
    if(visualize)
    {
        string PCname = "Plotting Points dx=" + to_string(g.dx()[0]);
        polyscope::registerPointCloud(PCname, s.xp());
        polyscope::getPointCloud(PCname)->addScalarQuantity("uplot", uplot);
        polyscope::getPointCloud(PCname)->addScalarQuantity("uexact", uexact);
        polyscope::getPointCloud(PCname)->addScalarQuantity("error", error);
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

    // ShiftedPoissonSphere(6, visualize);  // uncomment if you want to run just a single example with a single dx value 

    // make a lambda for convergence study, we need function format double f(double)
    auto Poisson = [&visualize](int& n)
    {
        return ShiftedPoissonSphere(n, visualize);
    };

    double avg_conv_order = h.ConvergenceStudy(4, 3, Poisson);
    if(abs(avg_conv_order - 2.1024515414334681) < 1e-16)
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
