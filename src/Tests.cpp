#include "Tests.h"

Tests::Tests()
{

}


// Test for standard interpolation matrix
void Tests::StandardInterpolationOrder()
{
    int numInterpDeg = 5;
    double avgOrder;
    for(int p = 0; p < numInterpDeg; ++p)
    {
        auto Interp2D = [this, &p](int& n)
        {
            return InterpolationOrder(n, p);
        };

        auto Interp3D = [this, &p](int& n)
        {
            return InterpolationOrder3D(n, p);
        };

        cout << "2D test with p = " << p << endl;
        h.ConvergenceStudy(6, 5, Interp2D);
        cout << " " << endl;
        cout << "3D test with p = " << p << endl;
        h.ConvergenceStudy(4, 3, Interp3D);
        cout << " " << endl;
    }
}

// Do a cp extension from a circle using interpolation and compare it to the exact cp extension
// Expected behaviour is to give p+1 order convergence, where p is the degree of the interpolating polynomial
vector<double> Tests::InterpolationOrder(int n, int interpPolyDegree)
{
    Grid g(interpPolyDegree, 2);
    Surface s("Circle");
    SpMat L, E;

    h.SetupMatrices(n, g, s, E, L);
    VectorXd u = h.SetInitialCondition(s);
    
    // apply cp extension
    VectorXd Eu = E * u;

    // compute error
    VectorXd error = Eu - u;

    vector<double> error_info(2);
    error_info[0] = error.lpNorm<Infinity>();
    error_info[1] = g.dx[0];
    return error_info;
}


vector<double> Tests::InterpolationOrder3D(int n, int interpPolyDegree)
{
    Grid g(interpPolyDegree, 3);
    Surface s("Sphere");
    SpMat L, E;
    
    h.SetupMatrices(n, g, s, E, L);
    VectorXd u = h.SetInitialCondition3D(s);

    // apply cp extension
    VectorXd Eu = E * u;

    // compute error
    VectorXd error = Eu - u;

    vector<double> error_info(2);
    error_info[0] = error.lpNorm<Infinity>();
    error_info[1] = g.dx[0];
    return error_info;
}