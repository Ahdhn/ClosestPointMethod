#pragma once

#include <vector>
#include <iostream>
#include "Node.h"

using namespace std;

// TO DO: Should use private? Everything in my classes are public, could run into problems where someone changes e.g., dx and it messes up the rest of the computation. Once it is set at the start none of this grid should be changed
class Grid {
    public:
        // Grid();
        Grid(int interpPolyDegree, int dimension);
        void InitializeGrid(int interpPolyDegree);

        // specify the bounding box
        vector<double> xstart; // (x, y, z,...) of negative corner of bounding box for rectangular grid
        vector<double> xend; // (x, y, z,...) of positive corner of bounding box for rectangular grid

        // grid spacing
        vector<double> dx; // dx, dy, dz,... to specify grid resolution. TO DO: current implementation only tested for dx=dy=dz

        int dim; // dimension of the embedding space
        int p; // interpolation order
        int order; // Laplacian order, used for bandwidth calculation. Might depend on gradient order or largest FD stencil

        vector<vector<double>> xg; // coordinates of the grid nodes
        vector<vector<double>> xb; // coordinates of the band nodes
        vector<Node> gridNodes; // grid nodes
        vector<Node> bandNodes; // grid nodes in the computational band

        vector<int> N1d; // number of grid points in each dimension
        int N; // total number of grid points

        // banding
        double bandwidth; // how many dx spacings should the bandwidth be
        int sizeBand; // total number of grid points in computational band

        void SetBoundingBox(vector<double> &xs, vector<double> &xe);
        void SetGridSpacing(vector<double> &deltax);

        void BuildRectangularGrid(int n);
        void SetNeighbours(Node &node, int index, int i, int j);
        void SetNeighbours(Node &node, int index, int i, int j, int k);

        void ComputeBandwidth();
        void SetComputationalBand(vector<double> &dist);

        void TestNeighbours(string gridName, bool isBanded);
        void TestNeighbours(string gridName);
        void TestNeighboursBand(string gridName);
};
