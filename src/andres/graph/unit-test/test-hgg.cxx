#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <vector>
#include <iterator>
#include <type_traits>

#include "andres/graph/graph.hxx"
#include "andres/graph/hyper-grid-graph.hxx"

inline void test(const bool& pred) {
    if (!pred) throw std::runtime_error("Test failed.");
}


// Define the OffsetVector type alias with a template parameter
template <unsigned char D>
using OffsetVector = std::vector<std::array<size_t, D>>;

//the HyperGridGraph type alias
template <unsigned char D>
using HyperGridGraph = andres::graph::HyperGridGraph<D>;

// Function to input and store the offset for a D-dimensional grid graph
template <unsigned char D>
OffsetVector<D> inputOffset() {
    OffsetVector<D> offset;
    size_t noffset;
    std::cout << "Enter the number of offset values for the D-dimensional grid graph:" << std::endl;
    std::cin >> noffset;

    for (size_t i = 0; i < noffset; ++i) {
        std::array<size_t, D> newArray;

        // Add elements to the array
        std::cout << "Enter the " << (i + 1) << " Offset list." << std::endl;
        for (size_t j = 0; j < D; ++j) {
            std::cout << "Enter the " << (j + 1) << " value:" << std::endl;
            size_t value;
            std::cin >> value;
            newArray[j] = value;
        }
        offset.push_back(newArray);
    }

    return offset;
}

// Function to print the offset values
template <unsigned char D>
void printOffset(const OffsetVector<D>& offset) {
    std::cout << "Offset values for the D-dimensional grid graph:" << std::endl;
    for (size_t i = 0; i < D; ++i) {
        for (const auto& value : offset[i]) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}

template <unsigned char D>
void test_empty_constructor() {
    // Create an empty hyper grid graph
    andres::graph::HyperGridGraph<D> hgg;

    // Check that the number of vertices is 0
    assert(hgg.numberOfVertices() == 0);

    std::cout << "Empty constructor test passed!" << std::endl;
}

void testEdge()

{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 0}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        graph({ 2, 4 },offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    
    std::cout << "Edge Coordinate test \n\n";
    for (size_t idx = 0; idx < graph.numberOfEdges(); ++idx)
    {

        graph.edge(idx, edge_coordinate);

        std::cout << "idx = " <<
            idx << ", pivot = (" <<

            edge_coordinate.pivotCoordinate_[0] <<
            ", " <<

            edge_coordinate.pivotCoordinate_[1] <<
            "), offset = " <<

            edge_coordinate.offsetIndex_ << "\n";

    }

}

void testVertex() //list of all vertex indices and their corresponding coordinates

{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 0}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        graph({ 2, 4 },offsetVector);
  
    
    andres::graph::HyperGridGraph<2>::VertexCoordinate 
        vertexC;


    std::cout << "Vertex Coordinate test \n\n";

    
    for (size_t idx = 0; idx < graph.numberOfVertices(); ++idx)
    {

        graph.vertex(idx, vertexC);

        std::cout << "idx = " <<
            idx << ", Vertex Coordinate = (" <<

            vertexC[0] <<
            ", " <<

            vertexC[1] <<")\n";

    }

}

int main()
{
   
    //const unsigned char D = 2; // Specify the dimensionality of the grid graph

    // Input the offset values based on the dimension D
    //OffsetVector<D> offset = inputOffset<D>();

    // Print the stored offset values
    //printOffset(offset);

    // Perform other tests
    testEdge();
    testVertex();
    //test_empty_constructor<D>();
    //test_empty_offsets();
    //test_vertex_of_edge();
    //test_find_edge();
    //test_iterate_vertices();
    //test_iterate_edges();
    //test_iterate_adjacent_vertices_of_vertex();
    // ...

    return 0;
}
