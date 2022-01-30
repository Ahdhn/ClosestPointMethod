#include "Grid.h"
#include "math.h"

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"


Grid::Grid(int interpPolyDegree, int dimension)
{
    InitializeGrid(interpPolyDegree);
    dim = dimension;
}


void Grid::InitializeGrid(int interpPolyDegree)
{
    p = interpPolyDegree;
    order = 2;
}


void Grid::SetBoundingBox(vector<double> &xs, vector<double> &xe)
{
    xstart = xs;
    xend = xe;
}


void Grid::SetGridSpacing(vector<double> &deltax)
{
    dx = deltax;
}


// build a rectangular grid
void Grid::BuildRectangularGrid(int n)
{
    // find number of grid points
    N1d.resize(dim);
    N = 1;
    for(int d = 0; d < dim; ++d)
    {
        N1d[d] = pow(2.0, n) + 1;
        N *= N1d[d];
    }

    // now build the uniform rectangular grid, with row-major lexicographical ordering (https://en.wikipedia.org/wiki/Row-_and_column-major_order)
    // TO DO: make this work for n dimension, see website above for indexing
    gridNodes.resize(N);
    xg.resize(N);
    int index;
    if(dim == 2)
    {
        for(int i = 0; i < N1d[0]; ++i)
        {
            for(int j = 0; j < N1d[1]; ++j)
            {
                index = i * N1d[1] + j; // shoots rays up +y direction, then moves over one x when a y line is done
                gridNodes[index].recIndex = index;
                xg[index].resize(dim);
                xg[index][0] = xstart[0] + i * dx[0];
                xg[index][1] = xstart[1] + j * dx[1];
            }
        }

        // set neighbours of the grid nodes
        for(int i = 0; i < N1d[0]; ++i)
        {
            for(int j = 0; j < N1d[1]; ++j)
            {
                index = i * N1d[1] + j;
                SetNeighbours(gridNodes[index], index, i, j);
            }
        }
    }
    else if(dim == 3)
    {
        for(int i = 0; i < N1d[0]; ++i)
        {
            for(int j = 0; j < N1d[1]; ++j)
            {
                for(int k = 0; k < N1d[2]; ++k)
                {
                    index = i * N1d[1] * N1d[2] + j * N1d[2] + k; // shoots rays in +z direction (from back to front), then moves up one y, then over one x when an yz-sheet is done
                    gridNodes[index].recIndex = index;
                    xg[index].resize(dim);
                    xg[index][0] = xstart[0] + i * dx[0];
                    xg[index][1] = xstart[1] + j * dx[1];
                    xg[index][2] = xstart[2] + k * dx[2];
                }
            }
        }

        // set neighbours of the grid nodes
        for(int i = 0; i < N1d[0]; ++i)
        {
            for(int j = 0; j < N1d[1]; ++j)
            {
                for(int k = 0; k < N1d[2]; ++k)
                {
                    index = i * N1d[1] * N1d[2] + j * N1d[2] + k; 
                    SetNeighbours(gridNodes[index], index, i, j, k);
                }
            }
        }
    }
}


// 2D version
void Grid::SetNeighbours(Node &node, int index, int i, int j)
{
    // right neighbour
    if(i + 1 < N1d[0])
    {
        node.neighbours[0] = &gridNodes[index + N1d[1]];
    } // nullptr otherwise set in Node constructor

    // left neighbour
    if(i - 1 >= 0)
    {
        node.neighbours[1] = &gridNodes[index - N1d[1]];
    }
    
    // up neighbour
    if( j + 1 < N1d[1])
    {
        node.neighbours[2] = &gridNodes[index + 1];
    }

    // down neighbour
    if(j - 1 >= 0)
    {
        node.neighbours[3] = &gridNodes[index - 1];
    }
}


// 3D version
void Grid::SetNeighbours(Node &node, int index, int i, int j, int k)
{
    // right neighbour
    if(i + 1 < N1d[0])
    {
        node.neighbours[0] = &gridNodes[index + N1d[1] * N1d[2]];
    }// nullptr otherwise set in Node constructor

    // left neighbour
    if(i - 1 >= 0)
    {
        node.neighbours[1] = &gridNodes[index - N1d[1] * N1d[2]];
    }
    
    // up neighbour
    if( j + 1 < N1d[1])
    {
        node.neighbours[2] = &gridNodes[index + N1d[2]];
    } 

    // down neighbour
    if(j - 1 >= 0)
    {
        node.neighbours[3] = &gridNodes[index - N1d[2]];
    }
    
    // front neighbour
    if( k + 1 < N1d[2])
    {
        node.neighbours[4] = &gridNodes[index + 1];
    }

    // back neighbour
    if(k - 1 >= 0)
    {
        node.neighbours[5] = &gridNodes[index - 1];
    }
}


//////////////////////////////////////////////////////////////////////////////
//    BANDING
//////////////////////////////////////////////////////////////////////////////


void Grid::ComputeBandwidth()
{
    bandwidth = 1.0001 * sqrt( (dim-1) * pow(0.5 * (p+1), 2) + pow(0.5 * (order + p + 1), 2) );
}


void Grid::SetComputationalBand(vector<double> &dist)
{
    if(dist.size() != N)
    {
        cout << "wrong size distance function" << endl;
    }

    sizeBand = 0;
    for(int i = 0; i < N; ++i)
    {
        if(abs(dist[i]) <= bandwidth * dx[0]) // TO DO: make work for different dx in each dimension
        {
            gridNodes[i].bandIndex = sizeBand;
            ++sizeBand;
            bandNodes.push_back(gridNodes[i]);
            xb.push_back(xg[i]);
        } // bandIndex = -1 otherwise, set in Node constructor
    }
    
    // remove neighbours not in band from bandNodes
    for(int i = 0; i < sizeBand; ++i)
    {
        for(int nbr = 0; nbr < 2 * dim; ++nbr)
        {
            if(bandNodes[i].neighbours[nbr]->bandIndex == -1)
            {
                bandNodes[i].neighbours[nbr] = nullptr;
            }
            else // set pointer to the bandNodes, instead of the gridNodes done during bandNodes.push_back(gridNodes[i]) above
            {
                bandNodes[i].neighbours[nbr] = &bandNodes[bandNodes[i].neighbours[nbr]->bandIndex];
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
//    TESTING
//////////////////////////////////////////////////////////////////////////////

// TO DO: Move testing to Tests.cpp

// visualize grid
void Grid::TestNeighbours(string gridName, bool isBanded)
{
    if(isBanded)
    {
        TestNeighboursBand(gridName);
    }
    else
    {
        TestNeighbours(gridName);
    }
}


// TO DO: visualization of this with polyscope
void Grid::TestNeighbours(string gridName)
{
    vector<int> startInd;
    vector<int> endInd;
    vector<double> edgeWeights;


    vector<double> xVertices(N);
    vector<double> yVertices(N);

    for(int i = 0; i < N; ++i)
    {
        xVertices[i] = xg[i][0];
        yVertices[i] = xg[i][1];

        for(int nbr = 0; nbr < 2 * dim; ++nbr)
        {
            if(gridNodes[i].neighbours[nbr] != nullptr)
            {
                edgeWeights.push_back(1.0);
                startInd.push_back(i + 1); // plus 1 for Matlab indexing
                endInd.push_back(gridNodes[i].neighbours[nbr]->recIndex + 1);
            }
        }
    }
}


void Grid::TestNeighboursBand(string gridName)
{
    vector<vector<double>> nodes(sizeBand, vector<double>(3));
    vector<double> pointCloudColor(sizeBand);
    vector<int> edge(2);
    vector<vector<int>> edges;
    for(int i = 0; i < sizeBand; ++i)
    {
        nodes[i][0] = xb[i][0];
        nodes[i][1] = xb[i][1];
        nodes[i][2] = 0;
        pointCloudColor[i] = nodes[i][2];

        for(int nbr = 0; nbr < 2 * dim; ++nbr)
        {
            if(bandNodes[i].neighbours[nbr] != nullptr)
            {
                edge[0] = i;
                edge[1] = bandNodes[i].neighbours[nbr]->bandIndex;
                edges.push_back(edge);
            }
        }
    }
    
    polyscope::registerPointCloud(gridName + "_PointCloud", nodes);
    polyscope::getPointCloud(gridName + "_PointCloud")->addScalarQuantity("whichGrid", pointCloudColor);
    
    polyscope::registerCurveNetwork(gridName + "_CurveNetwork", nodes, edges);
    polyscope::getCurveNetwork(gridName + "_CurveNetwork")->addNodeScalarQuantity("whichGrid", pointCloudColor);
}

