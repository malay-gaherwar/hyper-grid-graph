#pragma once
#ifndef ANDRES_GRAPH_GRID_GRAPH_HXX
#define ANDRES_GRAPH_GRID_GRAPH_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <iterator>
#include <array>
#include <initializer_list>
#include <vector>
#include "adjacency.hxx"
#include "visitor.hxx"

namespace andres {

namespace graph {



/// D-dimensional hyper grid graph.
template <unsigned char D = 2, class VISITOR = IdleGraphVisitor<std::size_t> >
class HyperGridGraph {
public:
    typedef std::size_t size_type;
    typedef VISITOR Visitor;
    typedef std::vector<std::array<size_type, D>> OffsetVector;      //definition of OffsetVector
    typedef andres::graph::Adjacency<size_type> AdjacencyType;

    static const size_type DIMENSION = static_cast<size_type> (D);

    typedef std::array<size_type, DIMENSION> VertexCoordinate;

    /// \class EdgeCoordinate
    /// Describes an edge as the integer index of the minimum of the
    /// two endpoints and the direction along which it is drawn.
    struct EdgeCoordinate {
        EdgeCoordinate(const VertexCoordinate&, const size_type, bool = false);
        EdgeCoordinate();
        /// The minimum of the two endpoints of the edge is specified
        /// as an integer and accessed by the \b pivot member.
        /// This specifies the \e origin of the edge, in a direction
        /// towards the positive orthant.
        /// (i.e.: in specifying the pivot of an edge, an isSmaller
        /// of false is implied.)
        VertexCoordinate pivotCoordinate_;
        /// The index of the offset along which the edge is drawn.
        size_type offset_index_;
    };

    HyperGridGraph(const Visitor & = Visitor());
    HyperGridGraph(const VertexCoordinate&, const OffsetVector&, const Visitor & = Visitor());
    HyperGridGraph(const std::initializer_list<std::size_t>, const Visitor & = Visitor());//not implemented yet
    void assign(const Visitor & = Visitor());
    void assign(const VertexCoordinate&, const Visitor & = Visitor());



    size_type numberOfEdgesFromVertex(const size_type) const;





private:
    // Member variables
    VertexCoordinate shape_;
    std::array<size_type, DIMENSION> edgeIndexOffsets_;
    std::array<size_type, DIMENSION> vertexIndexOffsets_;
    std::array<VertexCoordinate, DIMENSION> edgeShapes_;
    size_type numberOfVertices_;
    Visitor visitor_;

};



/// Construct an empty grid graph.
/// \tparam S the type of the indices used. It must be an unsigned type
/// (otherwise a compilation error is caused). Defaults to \c std::size_t .
/// \tparam VISITOR a visitor class to be used. Since this class is immutable
/// appart from resizing, this is a dummy variable meant for compatibility,
/// defaulting to \c andres::graph::IdgeGraphVisitor.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<unsigned char D, class VISITOR>
inline
    HyperGridGraph<D, VISITOR>::HyperGridGraph(
        const Visitor& visitor
    )
    : HyperGridGraph(VertexCoordinate({}), visitor) // Chain-call Constructor
{}


/// Construct a grid graph with a specified shape and offset.
///
/// \param shape the shape of the Grid Graph as an std::array.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<unsigned char D, class VISITOR>
inline
    HyperGridGraph<D, VISITOR>::HyperGridGraph(
        const VertexCoordinate& shape,
        const OffsetVector& offsets,
        const Visitor& visitor
    ) {
    assign(shape, offsets, visitor);
}






// Clear a grid graph. Assign function
///
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
///
template<unsigned char D, class VISITOR>
inline void
    HyperGridGraph<D, VISITOR>::assign(
        const Visitor& visitor
    ) {
    OffsetVector offset
    VertexCoordinate shape;
    std::fill(shape.begin(), shape.end(), 0);
    std::fill(offset.begin(), shape.end(), 0);
    assign(shape,offset, visitor);
}

/// Clear a grid graph and assign a new shape.
///
/// \param shape the shape of the grid graph
/// \param 
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
 /*
template<unsigned char D, class VISITOR> // doesnt work
inline void
HyperGridGraph<D, VISITOR>::assign(
    const VertexCoordinate& shape,
    const OffsetVector& offsets,
    const Visitor& visitor
) {
    shape_ = shape;
    offsets_ = offsets;
    visitor_ = visitor;

    // Set vertex offsets for fast vertex indexing.
    size_type cumprod = 1;
    for (size_type i = 0; i < D; ++i) {
        vertexIndexOffsets_[i] = cumprod;
        cumprod *= shape_[i];
    }
    numberOfVertices_ = cumprod; // Calculate the total number of vertices

    // Set edge offsets.
    size_type edgeIndexOffset = 0;
    for (size_type i = 0; i < D; ++i) {
        std::array<int, D>& edgeShape = edgeShapes_[i]; // Get the current edge shape array
        edgeShape = offsets_[i]; // Use the provided offsets for edge shapes

        // Adjust edge shapes by reducing each value by 1 if it's greater than 0
        for (size_type j = 0; j < D; ++j) {
            if (edgeShape[j] > 0) {
                --edgeShape[j];
            }
        }

        // Calculate the cumulative product of edgeShape values
        size_type cumprod = edgeShape[0];
        for (size_type j = 1; j < D; ++j) {
            cumprod *= edgeShape[j];
        }

        // Store the cumulative edge index offset and update edgeIndexOffset
        edgeIndexOffsets_[i] = (edgeIndexOffset += cumprod);
    }
}

 */       


//Calculation number of edges from vertex //completed
template<unsigned char D, class VISITOR> 
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::numberOfEdgesFromVertex(
    const size_type vertex
) const {
assert(vertex < numberOfVertices());
VertexCoordinate edgeCoordinate;
this->vertex(vertex, edgeCoordinate);
size_type numEdgesFromVertex = 0;

for (const auto& offset : offsets) {
    VertexCoordinate neighborCoordinate = edgeCoordinate;
    VertexCoordinate negNeighborCoordinate = edgeCoordinate;

    // Add the elements of the current offset array to the edge coordinate
    for (unsigned char i = 0; i < D; ++i) {
        neighborCoordinate[i] += offset[i];
        neighborCoordinateNeg[i] -= offset[i];
    }

    // Check if the neighbor coordinate is within the shape of the graph
    bool isValidNeighbor1 = true;
    bool isValidNeighbor2 = true;
    for (unsigned char i = 0; i < D; ++i) {
        if (neighborCoordinate[i] < 0 || neighborCoordinate[i] >= shape_[i]) {
            isValidNeighbor1 = false;
            //break;
            if (neighborCoordinateNeg[i] < 0 || neighborCoordinateNeg[i] >= shape_[i]) {
                isValidNeighbor2 = false;
                //break;
            }
        }

        if (isValidNeighbor1) {
            ++numEdgesFromVertex;
        }
        if (isValidNeighbor2) {
            ++numEdgesFromVertex;
        }

    }

    return numEdgesFromVertex;
}
/// Give Vertex Coordinate and retrieve the vertex index. //completed
/// \param[in] vertexIndex the integer index of the requested vertex
/// \param[out] vertexCoordinate The coordinates of the vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<unsigned char D, class VISITOR> 
void
GridGraph<D, VISITOR>::vertex(
    size_type vertexIndex,
    VertexCoordinate & vertexCoordinate
) const {
    assert(vertexIndex < numberOfVertices_);
    size_type i;
    for (i = 0; i < DIMENSION - 1; ++i) {
        vertexCoordinate[i] = vertexIndex % shape_[i];
        vertexIndex = vertexIndex / shape_[i];
    }
    vertexCoordinate[i] = vertexIndex;
}

/// Give Vertex index and retrieve the vertex coordinate. //complted
/// \param vertexCoordinate the integer index of the requested vertex
/// \return The integer index of the specified vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<unsigned char D, class VISITOR>
typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::vertex(
    const VertexCoordinate& vertexCoordinate
) const {
    size_type index = vertexCoordinate[DIMENSION - 1];
    for (size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * shape_[i - 1] + vertexCoordinate[i - 1];
    }
    return index;
}

/// Get the number of vertices.

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::numberOfVertices() const {
return numberOfVertices_;

    

 }

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_HYPER_GRID_GRAPH_HXX

