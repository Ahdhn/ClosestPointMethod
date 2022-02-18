#pragma once

class Node {
    public:
        Node();

        int recIndex() 
        { 
            return m_rec_index; 
        }; 
        
        int bandIndex() 
        { 
            return m_band_index; 
        }; 
        
        Node* neighbour(int nbr_index) 
        { 
            return m_neighbours[nbr_index]; 
        };

        void setRecIndex(int rec_index){ m_rec_index = rec_index; }; 
        void setBandIndex(int band_index){ m_band_index = band_index; }; 
        void setNeighbour(int nbr_index, Node* neighbour){ m_neighbours[nbr_index] = neighbour; };

    private:
        int m_rec_index; // index in original rectangular grid
        int m_band_index; // index in the computational band, instead of index in original rectangular grid. Set to -1 if grid point is not in band

        Node *m_neighbours[6]; // neighbours of current node, order is right(+x), left(-x), up(+y), down(-y), front(+z), back(-z). Set in Grid::SetNeighbours()
};
