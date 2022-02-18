#pragma once

#include <vector>
#include <iostream>
#include "Node.h"

using namespace std;

class Grid {
    public:
        Grid(int interp_deg, int embed_dim, double (&bBox)[2], int n);
        
        void SetComputationalBand(const vector<double> &dist);
        
        int dim(){ return m_embed_dim; };
        int p(){ return m_interp_deg; };
        int sizeBand(){ return m_size_band; };
        vector<double> dx(){ return m_dx; };
        vector<double> xstart(){ return m_xstart; };
        int N1d(int index){ return m_N1d[index]; };
        
        Node gridNode(int index){ return m_grid_nodes[index]; };
        Node bandNode(int index){ return m_band_nodes[index]; };

        const vector<vector<double>>& xg() const
        { 
            return m_grid_coords; 
        };
        
        const vector<vector<double>>& x() const 
        { 
            return m_band_coords; 
        };

    private:

        void ComputeBandwidth();
        void ConstructRectangularGrid(double (&bBox)[2], int n);
        void SetBoundingBox(vector<double> &xstart, vector<double> &xend);
        void BuildRectangularGrid(int n);
        void SetNeighbours(Node &node, int index, int i, int j);
        void SetNeighbours(Node &node, int index, int i, int j, int k);

        // specify the bounding box
        vector<double> m_xstart; // (x, y, z,...) of negative (bottom-left in 2D) corner of bounding box for rectangular grid
        vector<double> m_xend; // (x, y, z,...) of positive (top-right in 2D) corner of bounding box for rectangular grid

        // grid spacing
        vector<double> m_dx; // dx, dy, dz,... to specify grid resolution. TO DO: current implementation only tested for dx=dy=dz

        int m_embed_dim; // dimension of the embedding space
        int m_interp_deg; // polynomial interpolation order
        int m_order; // Laplacian order, used for bandwidth calculation. Might depend on gradient order or largest FD stencil

        vector<vector<double>> m_grid_coords; // coordinates of the grid nodes
        vector<vector<double>> m_band_coords; // coordinates of the band nodes
        vector<Node> m_grid_nodes; // grid nodes
        vector<Node> m_band_nodes; // grid nodes in the computational band

        vector<int> m_N1d; // number of grid points in each dimension
        int m_N; // total number of grid points

        // banding
        double m_bandwidth; // how many dx spacings should the bandwidth be
        int m_size_band; // total number of grid points in computational band
        
        void TestNeighbours(string gridName, bool isBanded);
        void TestNeighbours(string gridName);
        void TestNeighboursBand(string gridName);
};
