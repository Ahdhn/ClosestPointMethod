#pragma once

class Node {
    public:
        Node();

        int recIndex; // index in original rectangular grid
        int bandIndex; // index in the computational band, instead of index in original rectangular grid. Set to -1 if grid point is not in band

        Node *neighbours[6]; // neighbours of current node, order is right(+x), left(-x), up(+y), down(-y), front(+z), back(-z). Set in Grid::SetNeighbours()
};
