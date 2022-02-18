#include "Node.h"

Node::Node()
{
    setBandIndex(-1); // initially all grid nodes are considered not to be in the band. Computational band is built after cp and distance to surface is computed

    for(int i = 0; i < 6; ++i) // TO DO: Make this work for nD, right now I only have an array enough for 3D. Also, this is wasteful in 2D since there are always 2 nullptr neighbours
    {
        m_neighbours[i] = nullptr;
    }
}



