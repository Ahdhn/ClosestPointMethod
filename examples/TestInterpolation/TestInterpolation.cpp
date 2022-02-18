#include <vector>
#include <iostream>

#include "Grid.h"
#include "Surface.h"
#include "Interpolation.h"
#include "Helpers.h"

using namespace std;
using namespace Eigen;

Helpers h;

// Do a cp extension from a circle using interpolation and compare it to the exact cp extension
// Expected behaviour is to give p+1 order convergence, where p is the degree of the interpolating polynomial
vector<double> InterpolationOrder(int n, int interpPolyDegree)
{
    Surface s("Circle");
    Grid g(interpPolyDegree /* interpolation order */, 2 /* embedding space dimension */, s.bBox, n);
    SpMat L, E;

    h.SetupMatrices(g, s, E, L);
    VectorXd u = h.SetInitialCondition(s);
    
    // apply cp extension
    VectorXd Eu = E * u;

    // compute error
    VectorXd error = Eu - u;

    vector<double> error_info(2);
    error_info[0] = error.lpNorm<Infinity>();
    error_info[1] = g.dx()[0];
    return error_info;
}

vector<double> InterpolationOrder3D(int n, int interpPolyDegree)
{
    Surface s("Sphere");
    Grid g(interpPolyDegree/* interpolation order */, 3 /* embedding space dimension */, s.bBox, n);
    SpMat L, E;
    
    h.SetupMatrices(g, s, E, L);
    VectorXd u = h.SetInitialCondition3D(s);

    // apply cp extension
    VectorXd Eu = E * u;

    // compute error
    VectorXd error = Eu - u;

    vector<double> error_info(2);
    error_info[0] = error.lpNorm<Infinity>();
    error_info[1] = g.dx()[0];
    return error_info;
}

int main(int argc, char** argv)
{
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    int num_interp_deg = 5;
    vector<double> avg_order{0.9491689099929647, 1.9612014781538158, 2.9895758345733738, 3.9933848850038483, 4.9856689672597598};
    vector<double> avg_order3D{0.8623657582639215, 2.0033357901379887, 2.9085721368066584, 3.9279126234880501, 4.9144004098372678};

    bool all_passed = true;

    for(int p = 0; p < num_interp_deg; ++p)
    {
        auto Interp2D = [&p](int& n)
        {
            return InterpolationOrder(n, p);
        };

        auto Interp3D = [&p](int& n)
        {
            return InterpolationOrder3D(n, p);
        };

        double avg_conv_order = h.ConvergenceStudy(6, 5, Interp2D);
        if(abs(avg_conv_order - avg_order[p]) < 1e-16)
        {
            cout << "PASSED: 2D average converage order test with p = " << p << endl;
        }
        else
        {
            all_passed = false;
            cout << "FAILED: 2D average converage order test with p = " << p << endl;
        }

        double avg_conv_order3D = h.ConvergenceStudy(4, 3, Interp3D);
        if(abs(avg_conv_order3D - avg_order3D[p]) < 1e-16)
        {
            cout << "PASSED: 3D average converage order test with p = " << p << endl;
        }
        else
        {
            all_passed = false;
            cout << "FAILED: 3D average converage order test with p = " << p << endl;
        }
    }

    if(all_passed)
    {
        cout << "ALL PASSED :)" << endl;
    }
    else
    {
        cout << "AT LEAST ONE FAILED :(" << endl;
    }

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}
