#include "Grid.h"
#include "math.h"

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"


Grid::Grid(int interp_deg, int embed_dim, double (&bBox)[2], int n)
{
    m_interp_deg = interp_deg;
    m_embed_dim = embed_dim;
    m_order = 2;
    ComputeBandwidth();
    ConstructRectangularGrid(bBox, n);
}


void Grid::SetComputationalBand(const vector<double> &dist)
{
    if(dist.size() != m_N)
    {
        cout << "wrong size distance function" << endl;
    }

    m_size_band = 0;
    for(int i = 0; i < m_N; ++i)
    {
        if(abs(dist[i]) <= m_bandwidth * m_dx[0]) // TO DO: make work for different m_dx in each m_embed_dimension
        {
            m_grid_nodes[i].setBandIndex(m_size_band);
            ++m_size_band;
            m_band_nodes.push_back(m_grid_nodes[i]);
            m_band_coords.push_back(m_grid_coords[i]);
        } // bandIndex = -1 otherwise, set in Node constructor
    }
    
    // remove neighbours not in band from m_band_nodes
    for(int i = 0; i < m_size_band; ++i)
    {
        for(int nbr = 0; nbr < 2 * m_embed_dim; ++nbr)
        {
            if(m_band_nodes[i].neighbour(nbr)->bandIndex() == -1)
            {
                m_band_nodes[i].setNeighbour(nbr, nullptr);
            }
            else // set pointer to the m_band_nodes, instead of the m_grid_nodes done during m_band_nodes.push_back(m_grid_nodes[i]) above
            {
                m_band_nodes[i].setNeighbour(nbr, &m_band_nodes[m_band_nodes[i].neighbour(nbr)->bandIndex()]);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void Grid::ComputeBandwidth()
{
    m_bandwidth = 1.0001 * sqrt( (m_embed_dim-1) * pow(0.5 * (m_interp_deg+1), 2) + pow(0.5 * (m_order + m_interp_deg + 1), 2) );
}


// Set up rectangular embedding space grid within bounding box
void Grid::ConstructRectangularGrid(double (&bBox)[2], int n)
{
    int pad = 5;
    assert(pow(2.0, n) - 2.0 * pad > 0);
    double dx = (bBox[1] - bBox[0]) / (pow(2.0, n) - 2.0 * pad); // add extra 2*pad*dx width to the bounding box, while making sure the number of cells is 2^n
    vector<double> dx_vec(m_embed_dim, dx);
    m_dx = dx_vec;

    vector<double> xstart(m_embed_dim, bBox[0] - pad * dx); 
    vector<double> xend(m_embed_dim, bBox[1] + pad * dx);
    SetBoundingBox(xstart, xend); // x and y coordinates of bounding box

    BuildRectangularGrid(n);
}


void Grid::SetBoundingBox(vector<double> &xstart, vector<double> &xend)
{
    m_xstart = xstart;
    m_xend = xend;
}


// build a rectangular grid
void Grid::BuildRectangularGrid(int n)
{
    // find number of grid points
    m_N1d.resize(m_embed_dim);
    m_N = 1;
    for(int d = 0; d < m_embed_dim; ++d)
    {
        m_N1d[d] = pow(2.0, n) + 1;
        m_N *= m_N1d[d];
    }

    // now build the uniform rectangular grid, with row-major lexicographical ordering (https://en.wikipedia.org/wiki/Row-_and_column-major_order)
    // TO DO: make this work for n m_embed_dimension, see website above for indexing
    m_grid_nodes.resize(m_N);
    m_grid_coords.resize(m_N);
    int index;
    if(m_embed_dim == 2)
    {
        for(int i = 0; i < m_N1d[0]; ++i)
        {
            for(int j = 0; j < m_N1d[1]; ++j)
            {
                index = i * m_N1d[1] + j; // shoots rays up +y direction, then moves over one x when a y line is done
                m_grid_nodes[index].setRecIndex(index);
                m_grid_coords[index].resize(m_embed_dim);
                m_grid_coords[index][0] = m_xstart[0] + i * m_dx[0];
                m_grid_coords[index][1] = m_xstart[1] + j * m_dx[1];
            }
        }

        // set neighbours of the grid nodes
        for(int i = 0; i < m_N1d[0]; ++i)
        {
            for(int j = 0; j < m_N1d[1]; ++j)
            {
                index = i * m_N1d[1] + j;
                SetNeighbours(m_grid_nodes[index], index, i, j);
            }
        }
    }
    else if(m_embed_dim == 3)
    {
        for(int i = 0; i < m_N1d[0]; ++i)
        {
            for(int j = 0; j < m_N1d[1]; ++j)
            {
                for(int k = 0; k < m_N1d[2]; ++k)
                {
                    index = i * m_N1d[1] * m_N1d[2] + j * m_N1d[2] + k; // shoots rays in +z direction (from back to front), then moves up one y, then over one x when an yz-sheet is done
                    m_grid_nodes[index].setRecIndex(index);
                    m_grid_coords[index].resize(m_embed_dim);
                    m_grid_coords[index][0] = m_xstart[0] + i * m_dx[0];
                    m_grid_coords[index][1] = m_xstart[1] + j * m_dx[1];
                    m_grid_coords[index][2] = m_xstart[2] + k * m_dx[2];
                }
            }
        }

        // set neighbours of the grid nodes
        for(int i = 0; i < m_N1d[0]; ++i)
        {
            for(int j = 0; j < m_N1d[1]; ++j)
            {
                for(int k = 0; k < m_N1d[2]; ++k)
                {
                    index = i * m_N1d[1] * m_N1d[2] + j * m_N1d[2] + k; 
                    SetNeighbours(m_grid_nodes[index], index, i, j, k);
                }
            }
        }
    }
}


// 2D version
void Grid::SetNeighbours(Node &node, int index, int i, int j)
{
    // right neighbour
    if(i + 1 < m_N1d[0])
    {
        node.setNeighbour(0, &m_grid_nodes[index + m_N1d[1]]);
    } // nullptr otherwise set in Node constructor

    // left neighbour
    if(i - 1 >= 0)
    {
        node.setNeighbour(1, &m_grid_nodes[index - m_N1d[1]]);
    }
    
    // up neighbour
    if( j + 1 < m_N1d[1])
    {
        node.setNeighbour(2, &m_grid_nodes[index + 1]);
    }

    // down neighbour
    if(j - 1 >= 0)
    {
        node.setNeighbour(3, &m_grid_nodes[index - 1]);
    }
}


// 3D version
void Grid::SetNeighbours(Node &node, int index, int i, int j, int k)
{
    // right neighbour
    if(i + 1 < m_N1d[0])
    {
        node.setNeighbour(0, &m_grid_nodes[index + m_N1d[1] * m_N1d[2]]);
    }// nullptr otherwise set in Node constructor

    // left neighbour
    if(i - 1 >= 0)
    {
        node.setNeighbour(1, &m_grid_nodes[index - m_N1d[1] * m_N1d[2]]);
    }
    
    // up neighbour
    if( j + 1 < m_N1d[1])
    {
        node.setNeighbour(2, &m_grid_nodes[index + m_N1d[2]]);
    } 

    // down neighbour
    if(j - 1 >= 0)
    {
        node.setNeighbour(3, &m_grid_nodes[index - m_N1d[2]]);
    }
    
    // front neighbour
    if( k + 1 < m_N1d[2])
    {
        node.setNeighbour(4, &m_grid_nodes[index + 1]);
    }

    // back neighbour
    if(k - 1 >= 0)
    {
        node.setNeighbour(5, &m_grid_nodes[index - 1]);
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


    vector<double> xVertices(m_N);
    vector<double> yVertices(m_N);

    for(int i = 0; i < m_N; ++i)
    {
        xVertices[i] = m_grid_coords[i][0];
        yVertices[i] = m_grid_coords[i][1];

        for(int nbr = 0; nbr < 2 * m_embed_dim; ++nbr)
        {
            if(m_grid_nodes[i].neighbour(nbr) != nullptr)
            {
                edgeWeights.push_back(1.0);
                startInd.push_back(i + 1); // plus 1 for Matlab indexing
                endInd.push_back(m_grid_nodes[i].neighbour(nbr)->recIndex() + 1);
            }
        }
    }
}


void Grid::TestNeighboursBand(string gridName)
{
    vector<vector<double>> nodes(m_size_band, vector<double>(3));
    vector<double> pointCloudColor(m_size_band);
    vector<int> edge(2);
    vector<vector<int>> edges;
    for(int i = 0; i < m_size_band; ++i)
    {
        nodes[i][0] = m_band_coords[i][0];
        nodes[i][1] = m_band_coords[i][1];
        nodes[i][2] = 0;
        pointCloudColor[i] = nodes[i][2];

        for(int nbr = 0; nbr < 2 * m_embed_dim; ++nbr)
        {
            if(m_band_nodes[i].neighbour(nbr) != nullptr)
            {
                edge[0] = i;
                edge[1] = m_band_nodes[i].neighbour(nbr)->bandIndex();
                edges.push_back(edge);
            }
        }
    }
    
    polyscope::registerPointCloud(gridName + "_PointCloud", nodes);
    polyscope::getPointCloud(gridName + "_PointCloud")->addScalarQuantity("whichGrid", pointCloudColor);
    
    polyscope::registerCurveNetwork(gridName + "_CurveNetwork", nodes, edges);
    polyscope::getCurveNetwork(gridName + "_CurveNetwork")->addNodeScalarQuantity("whichGrid", pointCloudColor);
}

