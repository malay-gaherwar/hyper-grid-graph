#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <vector>
#include <iterator>
#include <type_traits>
#include <chrono>

#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/hyper-grid-graph.hxx"
#include "andres/graph/shortest-paths.hxx"

inline void test(const bool& pred) {
    if (!pred) throw std::runtime_error("Test failed.");
}



void testConstructor()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 0}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
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
    


    if (g.numberOfVertices() == 4 * 6 * 2 * 3)
    {
    std::cout << g.numberOfVertices() << "\n"<<4 * 6 * 2 * 3<< "\n";
    }
    
}
void testnumberOfEdges()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 2}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
    auto v = g.numberOfEdges();
    std::cout << "Number of edges: " << v << std::endl;

    
}
void testEdge1() //tests whether given edge index outputs correst edge coordinate for a 2D HGG

{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 0}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 50, 50 },offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    
    std::cout << "Edge Coordinate test \n\n";
    for (size_t idx = 0; idx < g.numberOfEdges(); ++idx)
    {

        g.edge(idx, edge_coordinate);

        std::cout << "idx = " <<
            idx << ", pivot = (" <<

            edge_coordinate.pivotCoordinate_[0] <<
            ", " <<

            edge_coordinate.pivotCoordinate_[1] <<
            "), offset index= " <<

            edge_coordinate.offsetIndex_ << "\n";

    }

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
        offsetVector{{ 1, 2}, { 2,10 }};
    andres::graph::HyperGridGraph<2>
        g({ 200, 30 }, offsetVector);
    
    auto v = g.vertexFromVertex(20, 1);
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
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        g({ 3, 3 }, offsetVector);
    andres::graph::HyperGridGraph<2>::AdjacencyIterator it(g, 4);
    /*
    typename andres::graph::HyperGridGraph<2>::AdjacencyIterator it;
   typename andres::graph::HyperGridGraph<2>::AdjacencyIterator end;
    
    // Iterate through the adjacent vertices
    for (; it != it.end(); ++it) {
        size_type adjacentVertex = *it;
        std::cout << "Vertex " << vertex << " has an adjacent neighbor: " << adjacentVertex << std::endl;
    */
}

void compareConstructor()
{

    auto start = std::chrono::high_resolution_clock::now();
    andres::graph::HyperGridGraph<4>::OffsetVector
        offsetVector{{ 1, 0, 0, 0}, { 0,1,0,0 }, { 0,0,1,0 }, { 0,0,0,1 }};
    andres::graph::HyperGridGraph<4>
        graph({ 2, 40000 ,3,4 }, offsetVector);
    andres::graph::HyperGridGraph<4>::EdgeCoordinate
        edge_coordinate;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time of HGG: " << duration.count() << " seconds." << std::endl;

    start = std::chrono::high_resolution_clock::now();
    andres::graph::GridGraph<4>
        g({ 2, 40000 ,3,4 });
    andres::graph::GridGraph<4>::EdgeCoordinate
        edge_coordinate1;
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Execution time of GG: " << duration.count() << " seconds." << std::endl;



}
template <size_t Dim>
void compareConstructionTime() {
    double duration1[99]; // Array to store HyperGridGraph execution times
    double duration2[99]; // Array to store GridGraph execution times

    for (size_t i = 2; i <= 100; ++i) {

        std::vector <std::array<int, Dim>> offset(Dim);
        for (auto& arr : offset) {
            arr.fill(0);
        }
        std::array<size_t, Dim> gridShape;
        for (size_t j = 0; j < Dim; ++j) {
            gridShape[j] = i;
            offset[j][j] = 1;
        }
        auto start = std::chrono::high_resolution_clock::now();
        andres::graph::HyperGridGraph<Dim>::OffsetVector offsetVector(offset);
        andres::graph::HyperGridGraph<Dim> graph(gridShape, offsetVector);
        andres::graph::HyperGridGraph<Dim>::EdgeCoordinate edge_coordinate;
        auto end = std::chrono::high_resolution_clock::now();
        duration1[i - 2] = std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count();

        start = std::chrono::high_resolution_clock::now();
        andres::graph::GridGraph<Dim> g(gridShape);
        andres::graph::GridGraph<Dim>::EdgeCoordinate edge_coordinate1;
        end = std::chrono::high_resolution_clock::now();
        duration2[i - 2] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    for (int i = 0; i < 99; ++i) {
        std::cout << "Shape (" << 2 + i << "," << 2 + i << "): HGG - " << duration1[i] << " nanoseconds, GG - " << duration2[i] << " nanoseconds" << std::endl;
    }
}

void testShortestPath()
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 2, 2}, { 0,2 }};
    andres::graph::HyperGridGraph<2>
        g({ 2, 3 }, offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;

    std::deque<std::size_t> path;

    bool found = andres::graph::spsp(g, 0, 100, path);
    std::cout << found << std::endl;
    for (const size_t& element : path) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    //test(found == true);
    //test(path.size() == 3);
    /*
    std::cout << "Edge Coordinate test \n\n";
    for (size_t idx = 0; idx < g.numberOfEdges(); ++idx)
    {

        g.edge(idx, edge_coordinate);

        std::cout << "idx = " <<
            idx << ", pivot = (" <<

            edge_coordinate.pivotCoordinate_[0] <<
            ", " <<

            edge_coordinate.pivotCoordinate_[1] <<
            "), offset index= " <<

            edge_coordinate.offsetIndex_ << "\n";

    }  */
}

void compareHGGvsG() // this code uses the shortest path algorithm to compare run time of HGG and G
{
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        hgg({ 100, 50 }, offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    andres::graph::HyperGridGraph<2>::VertexCoordinate vCoord;
    size_t n = 5000;
    andres::graph::Graph<> g(n);

 
    for (size_t i = 0; i < hgg.numberOfEdges(); i++) {
        hgg.edge(i, edge_coordinate);
        vCoord[0] = edge_coordinate.pivotCoordinate_[0] + offsetVector[edge_coordinate.offsetIndex_][0];
        vCoord[1] = edge_coordinate.pivotCoordinate_[1] + offsetVector[edge_coordinate.offsetIndex_][1];
       // std::cout << edge_coordinate.pivotCoordinate_[0] << ", " << edge_coordinate.pivotCoordinate_[1]
       //     << " offset index= " << edge_coordinate.offsetIndex_ << " destination coord: "<<vCoord[0]<<", "<<vCoord[1]<<"\n";
        size_t vIndex1 = hgg.vertex(edge_coordinate.pivotCoordinate_);
        size_t vIndex2 = hgg.vertex(vCoord);
        g.insertEdge(vIndex1, vIndex2);

               
    }
  
    std::deque<std::size_t> path;
    auto start = std::chrono::high_resolution_clock::now();
    bool found = andres::graph::spsp(g, 0, 2000, path);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time of graph: " << duration.count() << " seconds." << std::endl;

    std::deque<std::size_t> path1;
    start = std::chrono::high_resolution_clock::now();
    found = andres::graph::spsp(hgg, 0, 2000, path1);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Execution time of HGG: " << duration.count() << " seconds." << std::endl;
}

void compareHGGvsGWithWeights() {
    andres::graph::HyperGridGraph<2>::OffsetVector
        offsetVector{{ 1, 1}, { 0,1 }};
    andres::graph::HyperGridGraph<2>
        hgg({ 100, 50 }, offsetVector);
    andres::graph::HyperGridGraph<2>::EdgeCoordinate
        edge_coordinate;
    andres::graph::HyperGridGraph<2>::VertexCoordinate vCoord;
    size_t n = 5000;
    andres::graph::Graph<> g(n);

    std::vector<float> edgeWeights(hgg.numberOfEdges());

    for (size_t i = 0; i < hgg.numberOfEdges(); i++) {
        hgg.edge(i, edge_coordinate);
        vCoord[0] = edge_coordinate.pivotCoordinate_[0] + offsetVector[edge_coordinate.offsetIndex_][0];
        vCoord[1] = edge_coordinate.pivotCoordinate_[1] + offsetVector[edge_coordinate.offsetIndex_][1];
        // std::cout << edge_coordinate.pivotCoordinate_[0] << ", " << edge_coordinate.pivotCoordinate_[1]
        //     << " offset index= " << edge_coordinate.offsetIndex_ << " destination coord: "<<vCoord[0]<<", "<<vCoord[1]<<"\n";
        size_t vIndex1 = hgg.vertex(edge_coordinate.pivotCoordinate_);
        size_t vIndex2 = hgg.vertex(vCoord);
        g.insertEdge(vIndex1, vIndex2);

        edgeWeights[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);


    }
    float distance = 0;
    std::deque<std::size_t> path;
    auto start = std::chrono::high_resolution_clock::now();
    andres::graph::spsp(g, 0, 3000, edgeWeights.begin(), path, distance);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time of graph: " << duration.count() << " seconds. \nDistance: " << distance << std::endl;

    distance = 0;
    std::deque<std::size_t> path1;
    start = std::chrono::high_resolution_clock::now();
    andres::graph::spsp(hgg, 0, 3000, edgeWeights.begin(), path, distance);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Execution time of HGG: " << duration.count() << " seconds.\nDistance: " << distance << std::endl;
}

int main()
{
  
    // Perform tests
    
    //testConstructor(); //working
    //testnumberOfVertices(); //working
    //testnumberOfEdges(); //working
    //testEdge1(); //working
    //testEdge2(); //working
    //testVertex1(); //working
    //testVertex2(); //working
    //testShape(); //working 
    
    //testNumberOfEdgesFromVertex(); // working
    //testNumberOfEdgesToVertex();  //working same as numedgefromvertex
    
    //testvertexFromVertex(); //
    //testvertexToVertex(); //
    //testAdjacencyFromVertex(); //working
    //testFindEdge(); //working
    //testInsertedge();
    // 
    // 
    //testAdjacencyIterator(); //not working


    
    testShortestPath(); //

    //COMPARISION TESTS
    
    //compareConstructor();
    //compareConstructionTime<3>(); //completed
    //compareHGGvsG(); 
    //compareHGGvsGWithWeights();
 
    //test_empty_offsets();
    //test_vertex_of_edge();
    //test_find_edge();
    //test_iterate_vertices();
    //test_iterate_edges();
    //test_iterate_adjacent_vertices_of_vertex();
    // ...

    return 0;
}
