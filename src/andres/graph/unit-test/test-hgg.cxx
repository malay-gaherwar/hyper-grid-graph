#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <vector>
#include <iterator>
#include <type_traits>
#include <fstream>
#include <chrono>

#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/hyper-grid-graph.hxx"
#include "andres/graph/shortest-paths.hxx"


typedef std::array<double, 5> Runtime;

inline void test(const bool& pred) {
    if (!pred) throw std::runtime_error("Test failed.");
}


void testEmptyOffsetConstructor()
{
    andres::graph::HyperGridGraph<3>::OffsetVector
        offsetVector{};
    andres::graph::HyperGridGraph<3>
        g({ 2, 3,100 }, offsetVector);
     test(g.numberOfVertices() == 600);
     test(g.numberOfEdges()>10000000);

}
void testConstructor()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
   
    test(g.numberOfVertices() == 6);
    test(g.numberOfEdges() ==6);
}

void testnumberOfVertices()
{   
    andres::graph::HyperGridGraph<4>::OffsetVector
        offsetVector{{ 5, 1,2,0}, { 2,2,1,4 }};
    andres::graph::HyperGridGraph<4>
        g({ 4,6,2,3 }, offsetVector);
    andres::graph::HyperGridGraph<4>::EdgeCoordinate
        edge_coordinate;

    test(g.numberOfVertices() == 4 * 6 * 2 * 3);    
}
void testnumberOfEdges()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);   
    test(g.numberOfEdges() == 5);    
}
void testEdge1() //tests whether given edge index outputs correst edge coordinate for a 2D HGG

{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 0}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 3, 4 }, offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    g.edge(0, edge_coordinate);
    test(edge_coordinate.pivotCoordinate_[0] == 0);
    test(edge_coordinate.pivotCoordinate_[1] == 0);
    test(edge_coordinate.offsetIndex_ == 0);

    g.edge(4, edge_coordinate);
    test(edge_coordinate.pivotCoordinate_[0] == 0);
    test(edge_coordinate.pivotCoordinate_[1] == 2);
    test(edge_coordinate.offsetIndex_ == 0);

}
void testEdge2()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 4 }, offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    edge_coordinate.pivotCoordinate_ = { 0,1 };
    edge_coordinate.offsetIndex_ = 0;
    auto d = g.edge(edge_coordinate);
    std::cout << "Edge index" << d << std::endl;

}
void testVertex1() //list of all vertex indices and their corresponding coordinates

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
void testVertex2() //Given a vertex coordinate, outsputs the index
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 4 }, offsetVector);
    andres::graph::HyperGridGraph<2>::VertexCoordinate
        vertexC;
    vertexC = { 1,2};
    auto d = g.vertex(vertexC);
    std::cout << "Vertex index" << d << std::endl;
}
void testShape()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 5 }, offsetVector);
    auto d = g.shape(1);
    std::cout << "Shape: " << d << std::endl;

}

void testNumberOfEdgesFromVertex()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    auto num = g.numberOfVertices();
    std::array<size_t, 6> check;
    check = { 2, 1, 3, 3, 1, 2 };
    //std::cout << "Num of vert " <<  num <<  std::endl;
    for (size_t i = 0; i < g.numberOfVertices(); i++) {
        auto s = g.numberOfEdgesFromVertex(i);
        test(s == check[i]);

    }
}
void testNumberOfEdgesToVertex()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
    auto v = g.numberOfEdgesToVertex(0);
    auto z = g.numberOfEdgesFromVertex(1);
    std::cout << "Number of edges to vertex: " << v <<"   "<< z << std::endl;

}
void testvertexFromVertex() {

    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 3, 5 }, offsetVector);
    
    auto v = g.vertexFromVertex(12, 1);
    std::cout << "Vertex index at edge XXX from vertex index ZZZ: "  << v << std::endl;
}
void testvertexToVertex()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
    auto v = g.vertexToVertex(0, 0);
    std::cout << "Vertex number 0  to vertex index 1: " << v << std::endl;
}
void testAdjacencyFromVertex()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);

    size_t v;
    size_t z;
    auto a = g.adjacencyFromVertex(3, 0);
    std::cout << "Vertex index at edge XXX from vertex index ZZZ: " << a.vertex() <<", "<< a.edge() << std::endl;
}
void testFindEdge() {
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 3, 3 }, offsetVector);

    auto d = g.findEdge(3,7);
    std::cout << "Edge index btwn is " << d.first << ",  "  << d.second << std::endl;
}
void testInsertedge()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 3, 3 }, offsetVector);

    auto d = g.insertEdge(3, 7);
    std::cout << "Edge index btwn is " << d << std::endl;
}

void testAdjacencyIterator()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 5}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 3,400 }, offsetVector);
    size_t vertexIndex = 100;
    auto begin = g.adjacenciesFromVertexBegin(vertexIndex);
    auto end = g.adjacenciesFromVertexEnd(vertexIndex);

    // Iterate through the adjacencies and print information.
    std::cout << "Adjacencies from vertex " << vertexIndex << ":" << std::endl;
    for (auto it = begin; it != end; ++it) {
        const auto& adjacency = *it;
        std::cout << "Edge Index: " << adjacency.edge() << ", Vertex index: " << adjacency.vertex() << std::endl;
    }
}

int main()
{
  
    // Perform tests
    
    testEmptyOffsetConstructor();
    testConstructor(); 
    testnumberOfVertices(); 
    testnumberOfEdges(); 
    //testEdge1(); 
    //testEdge2(); 
    //testVertex1(); 
    //testVertex2();
    //testShape(); 
    
    //testNumberOfEdgesFromVertex(); 
    //testNumberOfEdgesToVertex();  
    
    //testvertexFromVertex(); 
    //testvertexToVertex(); 
    //testAdjacencyFromVertex(); 
    //testFindEdge(); //working
    //testInsertedge();

    //testAdjacencyIterator(); 
   
    //testShortestPath(); 


    return 0;
}
