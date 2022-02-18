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

        void GetClosestPoints(const vector<vector<double>> &x);
        void BandClosestPoints(Grid &g);
        void GetSurfaceParameterization(int &Np);

        vector<int> bdy(){ return m_bdy; };
        vector<double> surfaceParams(){ return m_surface_params; };

        const vector<vector<double>>& cpx() const
        { 
            return m_cpx; 
        };

        const vector<double>& dist() const
        { 
            return m_dist; 
        };

        const vector<vector<double>>& xp() const
        { 
            return m_xp; 
        };

        vector<double> thetap(){ return m_thetap; };
        vector<double> phip(){ return m_phip; };
        
        // TO DO: make this private
        double bBox[2];
    
    private:

        vector<double> m_surface_params;
        string m_surface;
        
        // closest point functions
        vector<vector<double>> m_cpxg; // closest points to nodes in the full rectangular grid
        vector<vector<double>> m_cpx; // closest points to nodes in the computational band
        vector<double> m_dist;
        vector<int> m_bdyg; // tag for which boundary the cp belongs to, 0 if not a boundary cp
        vector<int> m_bdy; // tag for which boundary the cp in band belongs to, 0 if not a boundary cp

        void cpCircle(const vector<vector<double>> &x);
        void cpSphere(const vector<vector<double>> &x);
        void cpArc(const vector<vector<double>> &x);
        void cpPoint2D(const vector<vector<double>> &x);

        // surface parametrization for visualization
        vector<vector<double>> m_xp;
        vector<double> m_thetap; // for circle and arc
        vector<double> m_phip; // for sphere

        void paramCircle(int &Np);
        void paramSphere(int &Np);
        void paramArc(int &Np);
};
