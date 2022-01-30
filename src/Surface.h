#pragma once

#include <iostream>
#include <vector>
#include <string>
#include "Grid.h"
#include "math.h"

using namespace std;

class Surface {
    public:
        Surface(string surf);
        Surface(string surf, vector<double> &params);

        vector<double> surfaceParams;
        string surface;
        double bBox[2];

        void GetClosestPoints(vector<vector<double>> &x);
        void BandClosestPoints(Grid &g);

        // closest point functions
        vector<vector<double>> cpxg; // closest points to nodes in the full rectangular grid
        vector<vector<double>> cpx; // closest points to nodes in the computational band
        vector<double> dist;
        vector<int> bdyg; // tag for which boundary the cp belongs to, 0 if not a boundary cp
        vector<int> bdy; // tag for which boundary the cp in band belongs to, 0 if not a boundary cp

        void cpCircle(vector<vector<double>> &x);
        void cpSphere(vector<vector<double>> &x);
        void cpArc(vector<vector<double>> &x);
        void cpPoint2D(vector<vector<double>> &x);

        // surface parametrization for visualization
        vector<vector<double>> xp;
        vector<double> thetap; // for circle and arc
        vector<double> phip; // for sphere
        void GetSurfaceParameterization(int &Np);

        void paramCircle(int &Np);
        void paramSphere(int &Np);
        void paramArc(int &Np);
};
