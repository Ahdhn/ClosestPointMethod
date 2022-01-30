#include <vector>
#include <iostream>

#include "Grid.h"
#include "Surface.h"
#include "FDMatrices.h"
#include "Interpolation.h"
#include "Helpers.h"

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"

using namespace std;
using namespace Eigen;

static bool g_run_sim = false;
static polyscope::PointCloud* g_plotting_points;

Helpers h;
Grid g(3 /* interpolation order */, 3 /* embedding space dimension */);
Surface s("Sphere");
SpMat L, E;
VectorXd u;

double powGridCell = 5;
double dt;
bool isImplicit = false;
    
SparseLU<SpMat> solver;
SpMat A; 

int Np = 10000;
SpMat plotE;


void runTimeStep(bool visualizeOutput = true)
{
    if(!isImplicit) // forward Euler time stepping
    {
        u = u + dt * (L * u);
        u = E * u; // closest point extension
    }
    else // implicit Euler time stepping
    {
        u = solver.solve(u);
        if(solver.info()!=Success) {
            cout << "solver failed" << endl;
        }
    }
    VectorXd uplot = plotE * u;

    if(visualizeOutput)
    {
        std::pair<double,double> bounds(-1.0,1.0);
        g_plotting_points->addScalarQuantity("uplot", uplot)->setMapRange(bounds);
    }
}

void callback()
{
    ImGui::PushItemWidth(100);

    if(ImGui::Button("Start time-stepping"))
    {
        g_run_sim = !g_run_sim;
    }

    if(g_run_sim)
    {
        runTimeStep();
    }
}


void initialize(bool visualizeOutput = true)
{
    h.SetupMatrices(powGridCell, g, s, E, L);

    u = h.SetInitialCondition(s);

    // time stepping heat flow
    dt = 0.2 * g.dx[0] * g.dx[0]; // gives second-order with both explicit and implicit Euler, as expected
    
    if(isImplicit) // for implicit Euler
    {
        A = h.ImplicitEulerMatrix(L, E, dt);
        solver.compute(A);
        if(solver.info()!=Success) {
            cout << "decomposition failed" << endl;
        }
    }

    plotE.resize(Np, g.sizeBand);
    h.GetPlotInterpMatrix(g, s, plotE, Np);
    VectorXd uplot = plotE * u;

    // Visualize with polyscope
    if(visualizeOutput)
    {   
        string cloudName = "Surface with dx = " + to_string(g.dx[0]);
        g_plotting_points = polyscope::registerPointCloud(cloudName, s.xp);
        g_plotting_points->addScalarQuantity("uplot", uplot)->setEnabled(true);
        g_plotting_points->setPointRadius(0.015);
    }
}


int main(int argc, char** argv)
{
    polyscope::init();

    initialize(true);

    polyscope::state::userCallback = callback;
    
    polyscope::show();
}
