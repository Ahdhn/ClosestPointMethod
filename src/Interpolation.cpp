#include "Interpolation.h"

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

// #define INTERPOLATION_TESTS

Interpolation::Interpolation(Grid &g, const vector<vector<double>> &xquery)
{
    xq = xquery;
    Nq = xq.size();

    // find base point for each query point (lower left corner of interpolation stencil (hypercube))
    Ibpt.resize(Nq);
    FindInterpolationBasePoint(g);

    w.resize(g.dim(), vector<vector<double>>(Nq, vector<double>(g.p()+1)));
    BuildInterpolationWeights(g, w);
}


void Interpolation::BuildInterpolationMatrix(Grid &g, SpMat &E)
{
    vector<T> coeffs; // list of non-zeros coefficients, (row, col, value) triplets to be used to construct sparse matrix afterward

    if(g.dim() == 2)
    {
        int xindex;
        int findex;
        for(int q = 0; q < Nq; ++q)
        {
            xindex = Ibpt[q];
            findex = Ibpt[q];
            for(int i = 0; i < g.p()+1; ++i)
            {
                for(int j = 0; j < g.p()+1; ++j)
                {
                    coeffs.push_back(T(q, g.bandNode(findex).bandIndex(), w[1][q][j] * w[0][q][i]));
                    findex = g.bandNode(findex).neighbour(2)->bandIndex(); // move one neighbour up
                }
                findex = g.bandNode(xindex).neighbour(0)->bandIndex(); // move one neighbour right, from the current x bottom base point
                xindex = findex;
            }
        }
    }
    else if(g.dim() == 3)
    {
        int xindex;
        int yindex;
        int zindex;
        for(int q = 0; q < Nq; ++q)
        {
            xindex = Ibpt[q];
            yindex = Ibpt[q];
            zindex = Ibpt[q];
            for(int i = 0; i < g.p()+1; ++i)
            {
                for(int j = 0; j < g.p()+1; ++j)
                {
                    for(int k = 0; k < g.p()+1; ++k)
                    {
                        coeffs.push_back(T(q, g.bandNode(zindex).bandIndex(), w[2][q][k] * w[1][q][j] * w[0][q][i]));
                        zindex = g.bandNode(zindex).neighbour(4)->bandIndex(); // move +z one neighbour
                    }
                    zindex = g.bandNode(yindex).neighbour(2)->bandIndex(); // move one neighbour up
                    yindex = zindex;
                }
                zindex = g.bandNode(xindex).neighbour(0)->bandIndex(); // move one neighbour right, from the current x bottom base point
                yindex = zindex;
                xindex = zindex;
            }
        }
    }
    else
    {
        cout << "Dimension of embedding space not implemented in Interpolation::BuildBandedGriddedInterpolationMatrix" << endl;
    }
    E.setFromTriplets(coeffs.begin(), coeffs.end());

#ifdef INTERPOLATION_TESTS
    TestRowSumOne(E);
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      ALIGNING STENCIL AROUND QUERY POINT
//////////////////////////////////////////////////////////////////////////////////////////////////////////


// Asymmetric picture
// p=0:   B====
// p=1:   B====x
// p=2:   B    x====x
// p=3:   B    x====x    x
// p=4:   B    x    x====x    x
// p=5:   B    x    x====x    x    x
// etc
// Symmetric picture
// p=0: ==B==
// p=1:   B====x
// p=2:   B  ==x==  x
// p=3:   B    x====x    x
// p=4:   B    x  ==x==  x    x
// p=5:   B    x    x====x    x    x
// etc
// TO DO: I use symmetric stencil below, test asymmetric also
void Interpolation::FindInterpolationBasePoint(Grid &g) // TO DO: put safety checks in here to ensure nodes in stencils exist in band
{
#ifdef INTERPOLATION_TESTS
    polyscope::init();
    vector<vector<double>> nodes(g.sizeBand, vector<double>(3));
    for(int i = 0; i < g.sizeBand; ++i)
    {
        nodes[i][0] = g.xb[i][0];
        nodes[i][1] = g.xb[i][1];
        nodes[i][2] = g.bandNodes[i].whichGrid;
    }
    polyscope::registerPointCloud("CompBand", nodes);
#endif
    
    vector<int> I1d(g.dim());
    for(int q = 0; q < Nq; ++q)
    {
        // TO DO: optimize this to get the linear index for any dimension, should also be done in Grid::BuildRectangularGrid
        if(g.p() % 2 == 0)
        {
            for(int d = 0; d < g.dim(); ++d)
            {
                // change below round to floor will give asymmetric interpolation stencil, but this breaks the stencils for even p in 3D, had to increase bandwidth for it to work
                I1d[d] = round((xq[q][d] - g.xstart()[d]) / g.dx()[d]); // index of nearest lower left (actually, it is lower left after the shifting below right?) grid point in the interpolation stencil
            }
        }
        else
        {
            for(int d = 0; d < g.dim(); ++d)
            {
                I1d[d] = floor((xq[q][d] - g.xstart()[d]) / g.dx()[d]); // index of nearest lower left (actually, it is lower left after the shifting below right?) grid point in the interpolation stencil
            }
        }
        
        Ibpt[q] = 0;
        for(int d = 0; d < g.dim(); ++d)
        {
            int bptMult = 1;
            for(int l = d+1; l < g.dim(); ++l)
            {
                bptMult *= g.N1d(l); // store just the linear index
            }
            Ibpt[q] += bptMult * I1d[d];
        }
        // above gives the nearest point in terms of the full rectangular grid, now change to banded grid
        Ibpt[q] = g.gridNode(Ibpt[q]).bandIndex();

        // move neighbours from the base point above to get to the correct lower left corner of the interpolation stencil, the above only gives a nearest grid point to the query point
        int numShifts = (g.p() % 2 == 0) ? g.p()/2 : (g.p()-1)/2;
        for(int i = 0; i < numShifts; ++i) 
        {
            for(int d = 0; d < g.dim(); ++d)
            {
                Ibpt[q] = g.bandNode(Ibpt[q]).neighbour(2 * d + 1)->bandIndex(); // move negative one neighbour, g.p/2 times, in each direction
            }
        }

#ifdef INTERPOLATION_TESTS
        VisualizeInterpolationStencil(q, g, Ibpt[q], xq[q], queryWhichGrid[q]);
#endif
    }
    
#ifdef INTERPOLATION_TESTS
    polyscope::show();
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
//    COMPUTING WEIGHTS AND INTERPOLATION MATRIX
/////////////////////////////////////////////////////////////////////////////////////////


void Interpolation::InterpolationWeights1D(vector<double> &x, vector<double> &w)
{
    int n = x.size();

    // see Berrut & Trefethen 2004, p. 504 for more information about loop below
    w[0] = 1.0;
    for(int j = 1; j < n; ++j)
    {
        for(int k = 0; k <= j-1; ++k)
        {
            w[k] = (x[k] - x[j]) * w[k];
        }
        
        w[j] = 1.0;
        for(int k = 0; k <= j-1; ++k)
        {
            w[j] = (x[j] - x[k]) * w[j];
        }
    }

    for(int j = 0; j < n; ++j)
    {
        w[j] = 1.0 / w[j];
        if(isinf(w[j])){cout << "Division by zero in Interpolation::InterpolationWeights1D" << endl;}
    }
}


void Interpolation::BuildInterpolationWeights1D(Grid &g, vector<double> &x, double &xq_subset, vector<double> &w)
{
    InterpolationWeights1D(x, w); // these are the weights for this set of grid points, does not depend on the query points (could maybe optimize by not recomputing if you get the same interpolation stencil for a different query point?)
    
    // check to see if any x[j] = xq
    bool isInterpPoint = false;
    int interpPoint;
    for(int j = 0; j < g.p()+1; ++j)
    {
        if(xq_subset == x[j])
        {
            isInterpPoint = true;
            interpPoint = j;
        }
    }

    if(isInterpPoint)
    {
        // set all weights to zero except the one for xq = x[j]
        for(int j = 0; j < g.p()+1; ++j)
        {
            w[j] = 0.0;
        }
        w[interpPoint] = 1.0;
    }
    else
    {
        // add dependence on the query point
        double wqSum = 0.0;
        for(int j = 0; j < g.p()+1; ++j)
        {
            w[j] /= xq_subset - x[j];
            wqSum += w[j];
        }

        for(int j = 0; j < g.p()+1; ++j)
        {
            w[j] /= wqSum;
        }
    }
}


void Interpolation::BuildInterpolationWeights(Grid &g, vector<vector<vector<double>>> &w)
{
    // compute interpolation weights
    int index;
    vector<double> x(g.p()+1);
    for(int d = 0; d < g.dim(); ++d) // for each direction x,y(,z)
    {
        for(int q = 0; q < Nq; ++q)
        {
            index = Ibpt[q];
            x[0] = g.x()[Ibpt[q]][d];
            for(int i = 1; i < g.p()+1; ++i)
            {
                x[i] = g.x()[g.bandNode(index).neighbour(2 * d)->bandIndex()][d];
                index = g.bandNode(index).neighbour(2 * d)->bandIndex(); // set index to your neighbour, to move one neighbour at a time
            }
            BuildInterpolationWeights1D(g, x, xq[q][d], w[d][q]);
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
//    TESTING
//////////////////////////////////////////////////////////////////////////////

int intersection(std::vector<int> &v1, int &v2num)
{
    int index = -1;
    for(int i = 0; i < v1.size(); ++i)
    {
        if(v1[i] == v2num)
        {
            index = i;
            break;
        }
    }

    return index;
}


void Interpolation::VisualizeInterpolationStencil(int &q, Grid &g, int &Ibpt, vector<double> &xq, int &queryWG)
{
    int xindex;
    int yindex;
    int findex;

    vector<int> stencilIndices;
    xindex = Ibpt;
    findex = Ibpt;
    for(int i = 0; i < g.p()+1; ++i)
    {
        for(int j = 0; j < g.p()+1; ++j)
        {
            stencilIndices.push_back(g.bandNode(findex).bandIndex());
            if(g.bandNode(findex).neighbour(2) == nullptr)
            {
                cout << "problem with up neighbour in Interpolation::TestInterpolationStencil" << endl;
            }
            findex = g.bandNode(findex).neighbour(2)->bandIndex(); // move one neighbour up
        }
        if(g.bandNode(xindex).neighbour(0) == nullptr)
        {
            cout << "problem with right neighbour in Interpolation::TestInterpolationStencil" << endl;
        }
        findex = g.bandNode(xindex).neighbour(0)->bandIndex(); // move one neighbour right, from the current x bottom base point
        xindex = findex;
    }

    int numStencil = (g.p() + 1) * (g.p() + 1);
    vector<vector<double>> nodes(numStencil, vector<double>(3));
    int index;
    for(int i = 0; i < numStencil; ++i)
    {
        index = stencilIndices[i];
        nodes[i][0] = g.x()[index][0];
        nodes[i][1] = g.x()[index][1];
    }

    int matchingIndex;
    int endIndex;
    vector<int> edge(2);
    vector<vector<int>> edges;
    for(int i = 0; i < numStencil; ++i)
    {
        index = stencilIndices[i];
        for(int nbr = 0; nbr < 2 * g.dim(); ++nbr)
        {
            if(g.bandNode(index).neighbour(nbr) != nullptr)
            {
                endIndex = g.bandNode(index).neighbour(nbr)->bandIndex();
                matchingIndex = intersection(stencilIndices, endIndex);
                
                if(matchingIndex != -1)
                {
                    edge[0] = i;
                    edge[1] = matchingIndex;
                    edges.push_back(edge);
                }
            }
        }
    }

    vector<vector<double>> xquery(1, vector<double>(3));
    xquery[0][0] = xq[0];
    xquery[0][1] = xq[1];
    xquery[0][2] = queryWG;
    polyscope::registerPointCloud("QueryPoint" + to_string(q), xquery);
    polyscope::getPointCloud("QueryPoint" + to_string(q))->setEnabled(false);
    polyscope::registerCurveNetwork("CurveNetwork" + to_string(q), nodes, edges);
    polyscope::getCurveNetwork("CurveNetwork" + to_string(q))->setEnabled(false);
}


void Interpolation::TestRowSumOne(SpMat &E)
{
    Eigen::VectorXd sumVec = E * Eigen::VectorXd::Ones(E.cols());
    if(!sumVec.isOnes())
    {
        cout << "Interpolation matrix rows do not all sum to one" << endl;
    }
}
