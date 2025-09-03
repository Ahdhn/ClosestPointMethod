#define _USE_MATH_DEFINES 
#include <cmath>

#include "Surface.h"


// surfaces with default parameters
Surface::Surface(string surf)
{
    m_surface = surf;
    
    if(m_surface.compare("Circle") == 0)
    {
        m_surface_params.push_back(1.0); // radius of the circle
        bBox[0] = -m_surface_params[0];
        bBox[1] = m_surface_params[0];
    }

    if(m_surface.compare("Sphere") == 0)
    {
        m_surface_params.push_back(1.0); // radius of the sphere
        bBox[0] = -m_surface_params[0];
        bBox[1] = m_surface_params[0];
    }

    if(m_surface.compare("Arc") == 0)
    {
        m_surface_params.push_back(1.0); // radius of the circle
        m_surface_params.push_back(-3*M_PI_4); // angle1
        m_surface_params.push_back(M_PI_4); // angle2
        bBox[0] = -m_surface_params[0];
        bBox[1] = m_surface_params[0];
    }
}


// allow user to specify surface parameters instead
Surface::Surface(string surf, vector<double> &params)
{
    m_surface = surf;
    m_surface_params = params;
    bBox[0] = -m_surface_params[0];
    bBox[1] = m_surface_params[0];
}


void Surface::GetClosestPoints(const vector<vector<double>> &x)
{
    m_cpxg.resize(x.size(), vector<double>(x[0].size()));
    m_dist.resize(x.size());
    m_bdyg.resize(x.size(), 0.0);

    if(m_surface.compare("Circle") == 0)
    {
        cpCircle(x);
    }
    else if(m_surface.compare("Point2D") == 0)
    {
        cpPoint2D(x);
    }
    else if(m_surface.compare("Arc") == 0)
    {
        cpArc(x);
    }
    else if(m_surface.compare("Sphere") == 0)
    {
        cpSphere(x);
    }
    else
    {
        cout << "Invalid surface type in Surface::GetClosestPoints." << endl;
    }
}


// make a vector of just the closest points in the computational band
void Surface::BandClosestPoints(Grid &g)
{
    m_cpx.resize(g.sizeBand(), vector<double>(g.dim()));
    for(int i = 0; i < g.sizeBand(); ++i)
    {
        for(int d = 0; d < g.dim(); ++d)
        {
            m_cpx[i][d] = m_cpxg[g.bandNode(i).recIndex()][d];
        }
    }
    
    // for open surfaces, all zeros if it is a closed surface
    m_bdy.resize(g.sizeBand());
    for(int i = 0; i < g.sizeBand(); ++i)
    {
        m_bdy[i] = m_bdyg[g.bandNode(i).recIndex()];
    }
}

void Surface::GetSurfaceParameterization(int &Np)
{
    if(m_surface.compare("Circle") == 0)
    {
        m_xp.resize(Np, vector<double>(2));
        m_thetap.resize(Np);
        paramCircle(Np);
    }
    else if(m_surface.compare("Sphere") == 0)
    {
        m_xp.resize(Np, vector<double>(3));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramSphere(Np);
    }
    else if(m_surface.compare("Arc") == 0)
    {
        m_xp.resize(Np, vector<double>(2));
        m_thetap.resize(Np);
        paramArc(Np);
    }
    else
    {
        cout << "Invalid surface type in Surface::GetSurfaceParameterization." << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//          Private Functions
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
//          CLOSEST POINT FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////

// Closest point to a circle embedded in 2 or 3 dimensions
// TO DO: do it for n-dimensions
void Surface::cpCircle(const vector<vector<double>> &x)
{
    double R = m_surface_params[0];
    double theta;
    double gridNorm;
    int dim = x[0].size();
    for(int i = 0; i < x.size(); ++i)
    {
        theta = atan2(x[i][1], x[i][0]);
        m_cpxg[i][0] = R * cos(theta);
        m_cpxg[i][1] = R * sin(theta);
        if(dim == 3)
        {
            m_cpxg[i][2] = 0.0;
        }
        
        m_dist[i] = 0.0;
        for(int d = 0; d < dim; ++d)
        {
            m_dist[i] += pow(x[i][d] - m_cpxg[i][d], 2);
        }
        m_dist[i] = sqrt(m_dist[i]);
    }
}


void Surface::cpSphere(const vector<vector<double>> &x)
{
    double R = m_surface_params[0];
    double theta;
    double phi;

    for(int i = 0; i < x.size(); ++i)
    {
        theta = atan2(x[i][1], x[i][0]);
        m_dist[i] = sqrt( pow(x[i][0],2) + pow(x[i][1],2) + pow(x[i][2],2) );
        if(m_dist[i] != 0)
        {
            phi = acos(x[i][2] / m_dist[i]);
        }
        else
        {
            phi = 0.0; // pick a point on the surface for the (0,0,0) grid point
        }

        m_cpxg[i][0] = R * sin(phi) * cos(theta);
        m_cpxg[i][1] = R * sin(phi) * sin(theta);
        m_cpxg[i][2] = R * cos(phi);
        
        m_dist[i] -= R;
    }
}


// arc specified by portion of the circle between [angle1, angle2]. These angles must be specified in the interval (-pi, pi].
void Surface::cpArc(const vector<vector<double>> &x)
{
    double R = m_surface_params[0];
    double angle1 = m_surface_params[1];
    double angle2 = m_surface_params[2];
    int dim = x[0].size();

    if(angle1 < -M_PI)
    {
        cout << "angle1 in Surface::cpArc should be in (-pi, pi]." << endl;
    }

    if(angle2 > M_PI)
    {
        cout << "angle2 in Surface::cpArc should be in (-pi, pi]." << endl;
    }

    if(angle1 >= angle2)
    {
        cout << "angle1 should be less than angle2 in Surface::cpArc." << endl;
    }


    // find the point opposite the half point in the arc to split the boundary points
    double halfDist = 0.5 * (2 * M_PI - (angle2 - angle1));

    double theta;
    m_bdyg.resize(x.size());
    for(int i = 0; i < x.size(); ++i)
    {
        theta = atan2(x[i][1], x[i][0]);

        if(theta > angle1 & theta < angle2)
        {
            m_cpxg[i][0] = R * cos(theta);
            m_cpxg[i][1] = R * sin(theta);
            m_bdyg[i] = 0;
        }
        else 
        {
            if(theta <= angle1 & theta >= angle1 - halfDist)
            {
                m_bdyg[i] = 1; // assign point to first boundary
                m_cpxg[i][0] = R * cos(angle1);
                m_cpxg[i][1] = R * sin(angle1);
            }
            else if(theta <= angle1 & theta < angle1 - halfDist)
            {
                m_bdyg[i] = 2; // assign point to second boundary
                m_cpxg[i][0] = R * cos(angle2);
                m_cpxg[i][1] = R * sin(angle2);
            }
            else if(theta >= angle2 & theta < angle2 + halfDist)
            {
                m_bdyg[i] = 2; // assign point to second boundary
                m_cpxg[i][0] = R * cos(angle2);
                m_cpxg[i][1] = R * sin(angle2);
            }
            else if(theta >= angle2 & theta >= angle2 + halfDist)
            {
                m_bdyg[i] = 1; // assign point to second boundary
                m_cpxg[i][0] = R * cos(angle1);
                m_cpxg[i][1] = R * sin(angle1);
            }
        }

        if(dim == 3)
        {
            m_cpxg[i][2] = 0.0;
        }
        
        m_dist[i] = 0.0;
        for(int d = 0; d < dim; ++d)
        {
            m_dist[i] += pow(x[i][d] - m_cpxg[i][d], 2);
        }
        m_dist[i] = sqrt(m_dist[i]);
    }
}


void Surface::cpPoint2D(const vector<vector<double>> &x)
{
    double xc = m_surface_params[0];
    double yc = m_surface_params[1];

    for(int i = 0; i < x.size(); ++i)
    {
        m_cpxg[i][0] = xc;
        m_cpxg[i][1] = yc;
        
        m_dist[i] = sqrt(pow(x[i][0] - xc, 2) + pow(x[i][1] - yc, 2));
    }
}


///////////////////////////////////////////////////////////////////////////////////
//  SURFACE PARAMETRIZATIONS (for visualization - final solution on nicely spaced points)
///////////////////////////////////////////////////////////////////////////////////

void Surface::paramCircle(int &Np)
{
    double R = m_surface_params[0];
    double dtheta = 2 * M_PI / (Np - 1);

    for(int i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * dtheta;
        m_xp[i][0] = R * cos(m_thetap[i]);
        m_xp[i][1] = R * sin(m_thetap[i]);
    }
}


// Note: the number of parametrization points, Np, must be the square of some number here
void Surface::paramSphere(int &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        cout << "Np must be the square of some number for Surface::paramSphere" << endl;
    }
    int sqrtNp = sqrt(Np);
    double R = m_surface_params[0];
    double dtheta = 2 * M_PI / (sqrtNp - 1);
    double dphi = M_PI / (sqrtNp - 1);

    for(int i = 0; i < sqrtNp; ++i)
    {
        for(int j = 0; j < sqrtNp; ++j)
        {
            m_thetap[i * sqrtNp + j] = i * dtheta;
            m_phip[i * sqrtNp + j] = j * dphi;
        }
    }

    for(int i = 0; i < Np; ++i)
    {
        m_xp[i][0] = R * sin(m_phip[i]) * cos(m_thetap[i]);
        m_xp[i][1] = R * sin(m_phip[i]) * sin(m_thetap[i]);
        m_xp[i][2] = R * cos(m_phip[i]);
    }
}


void Surface::paramArc(int &Np)
{
    double R = m_surface_params[0];
    double angle1 = m_surface_params[1];
    double angle2 = m_surface_params[2];
    double dtheta = (angle2 - angle1) / (Np - 1);

    for(int i = 0; i < Np; ++i)
    {
        m_thetap[i] = angle1 + i * dtheta;
        m_xp[i][0] = R * cos(m_thetap[i]);
        m_xp[i][1] = R * sin(m_thetap[i]);
    }
}
