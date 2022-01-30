#include "Surface.h"


// surfaces with default parameters
Surface::Surface(string surf)
{
    surface = surf;
    
    if(surface.compare("Circle") == 0)
    {
        surfaceParams.push_back(1.0); // radius of the circle
        bBox[0] = -surfaceParams[0];
        bBox[1] = surfaceParams[0];
    }

    if(surface.compare("Sphere") == 0)
    {
        surfaceParams.push_back(1.0); // radius of the sphere
        bBox[0] = -surfaceParams[0];
        bBox[1] = surfaceParams[0];
    }

    if(surface.compare("Arc") == 0)
    {
        surfaceParams.push_back(1.0); // radius of the circle
        surfaceParams.push_back(-3*M_PI_4); // angle1
        surfaceParams.push_back(M_PI_4); // angle2
        bBox[0] = -surfaceParams[0];
        bBox[1] = surfaceParams[0];
    }
}


// allow user to specify surface parameters instead
Surface::Surface(string surf, vector<double> &params)
{
    surface = surf;
    surfaceParams = params;
    bBox[0] = -surfaceParams[0];
    bBox[1] = surfaceParams[0];
}


void Surface::GetClosestPoints(vector<vector<double>> &x)
{
    cpxg.resize(x.size(), vector<double>(x[0].size()));
    dist.resize(x.size());
    bdyg.resize(x.size(), 0.0);

    if(surface.compare("Circle") == 0)
    {
        cpCircle(x);
    }
    else if(surface.compare("Point2D") == 0)
    {
        cpPoint2D(x);
    }
    else if(surface.compare("Arc") == 0)
    {
        cpArc(x);
    }
    else if(surface.compare("Sphere") == 0)
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
    cpx.resize(g.sizeBand, vector<double>(g.dim));
    for(int i = 0; i < g.sizeBand; ++i)
    {
        for(int d = 0; d < g.dim; ++d)
        {
            cpx[i][d] = cpxg[g.bandNodes[i].recIndex][d];
        }
    }
    
    // for open surfaces, all zeros if it is a closed surface
    bdy.resize(g.sizeBand);
    for(int i = 0; i < g.sizeBand; ++i)
    {
        bdy[i] = bdyg[g.bandNodes[i].recIndex];
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//          CLOSEST POINT FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////


// Closest point to a circle embedded in 2 or 3 dimensions
// TO DO: do it for n-dimensions
void Surface::cpCircle(vector<vector<double>> &x)
{
    double R = surfaceParams[0];
    double theta;
    double gridNorm;
    int dim = x[0].size();
    for(int i = 0; i < x.size(); ++i)
    {
        theta = atan2(x[i][1], x[i][0]);
        cpxg[i][0] = R * cos(theta);
        cpxg[i][1] = R * sin(theta);
        if(dim == 3)
        {
            cpxg[i][2] = 0.0;
        }
        
        dist[i] = 0.0;
        for(int d = 0; d < dim; ++d)
        {
            dist[i] += pow(x[i][d] - cpxg[i][d], 2);
        }
        dist[i] = sqrt(dist[i]);
    }
}


void Surface::cpSphere(vector<vector<double>> &x)
{
    double R = surfaceParams[0];
    double theta;
    double phi;

    for(int i = 0; i < x.size(); ++i)
    {
        theta = atan2(x[i][1], x[i][0]);
        dist[i] = sqrt( pow(x[i][0],2) + pow(x[i][1],2) + pow(x[i][2],2) );
        if(dist[i] != 0)
        {
            phi = acos(x[i][2] / dist[i]);
        }
        else
        {
            phi = 0.0; // pick a point on the surface for the (0,0,0) grid point
        }

        cpxg[i][0] = R * sin(phi) * cos(theta);
        cpxg[i][1] = R * sin(phi) * sin(theta);
        cpxg[i][2] = R * cos(phi);
        
        dist[i] -= R;
    }
}


// arc specified by portion of the circle between [angle1, angle2]. These angles must be specified in the interval (-pi, pi].
void Surface::cpArc(vector<vector<double>> &x)
{
    double R = surfaceParams[0];
    double angle1 = surfaceParams[1];
    double angle2 = surfaceParams[2];
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
    bdyg.resize(x.size());
    for(int i = 0; i < x.size(); ++i)
    {
        theta = atan2(x[i][1], x[i][0]);

        if(theta > angle1 & theta < angle2)
        {
            cpxg[i][0] = R * cos(theta);
            cpxg[i][1] = R * sin(theta);
            bdyg[i] = 0;
        }
        else 
        {
            if(theta <= angle1 & theta >= angle1 - halfDist)
            {
                bdyg[i] = 1; // assign point to first boundary
                cpxg[i][0] = R * cos(angle1);
                cpxg[i][1] = R * sin(angle1);
            }
            else if(theta <= angle1 & theta < angle1 - halfDist)
            {
                bdyg[i] = 2; // assign point to second boundary
                cpxg[i][0] = R * cos(angle2);
                cpxg[i][1] = R * sin(angle2);
            }
            else if(theta >= angle2 & theta < angle2 + halfDist)
            {
                bdyg[i] = 2; // assign point to second boundary
                cpxg[i][0] = R * cos(angle2);
                cpxg[i][1] = R * sin(angle2);
            }
            else if(theta >= angle2 & theta >= angle2 + halfDist)
            {
                bdyg[i] = 1; // assign point to second boundary
                cpxg[i][0] = R * cos(angle1);
                cpxg[i][1] = R * sin(angle1);
            }
        }

        if(dim == 3)
        {
            cpxg[i][2] = 0.0;
        }
        
        dist[i] = 0.0;
        for(int d = 0; d < dim; ++d)
        {
            dist[i] += pow(x[i][d] - cpxg[i][d], 2);
        }
        dist[i] = sqrt(dist[i]);
    }
}


void Surface::cpPoint2D(vector<vector<double>> &x)
{
    double xc = surfaceParams[0];
    double yc = surfaceParams[1];

    for(int i = 0; i < x.size(); ++i)
    {
        cpxg[i][0] = xc;
        cpxg[i][1] = yc;
        
        dist[i] = sqrt(pow(x[i][0] - xc, 2) + pow(x[i][1] - yc, 2));
    }
}


///////////////////////////////////////////////////////////////////////////////////
//  SURFACE PARAMETRIZATIONS (for visualization - final solution on nicely spaced points)
///////////////////////////////////////////////////////////////////////////////////


void Surface::GetSurfaceParameterization(int &Np)
{
    if(surface.compare("Circle") == 0)
    {
        xp.resize(Np, vector<double>(2));
        thetap.resize(Np);
        paramCircle(Np);
    }
    else if(surface.compare("Sphere") == 0)
    {
        xp.resize(Np, vector<double>(3));
        thetap.resize(Np);
        phip.resize(Np);
        paramSphere(Np);
    }
    else if(surface.compare("Arc") == 0)
    {
        xp.resize(Np, vector<double>(2));
        thetap.resize(Np);
        paramArc(Np);
    }
    else
    {
        cout << "Invalid surface type in Surface::GetSurfaceParameterization." << endl;
    }
}


void Surface::paramCircle(int &Np)
{
    double R = surfaceParams[0];
    double dtheta = 2 * M_PI / (Np - 1);

    for(int i = 0; i < Np; ++i)
    {
        thetap[i] = i * dtheta;
        xp[i][0] = R * cos(thetap[i]);
        xp[i][1] = R * sin(thetap[i]);
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
    double R = surfaceParams[0];
    double dtheta = 2 * M_PI / (sqrtNp - 1);
    double dphi = M_PI / (sqrtNp - 1);

    for(int i = 0; i < sqrtNp; ++i)
    {
        for(int j = 0; j < sqrtNp; ++j)
        {
            thetap[i * sqrtNp + j] = i * dtheta;
            phip[i * sqrtNp + j] = j * dphi;
        }
    }

    for(int i = 0; i < Np; ++i)
    {
        xp[i][0] = R * sin(phip[i]) * cos(thetap[i]);
        xp[i][1] = R * sin(phip[i]) * sin(thetap[i]);
        xp[i][2] = R * cos(phip[i]);
    }
}


void Surface::paramArc(int &Np)
{
    double R = surfaceParams[0];
    double angle1 = surfaceParams[1];
    double angle2 = surfaceParams[2];
    double dtheta = (angle2 - angle1) / (Np - 1);

    for(int i = 0; i < Np; ++i)
    {
        thetap[i] = angle1 + i * dtheta;
        xp[i][0] = R * cos(thetap[i]);
        xp[i][1] = R * sin(thetap[i]);
    }
}
