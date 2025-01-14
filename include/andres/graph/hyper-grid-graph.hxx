#pragma once
#ifndef ANDRES_GRAPH_HYPER_GRID_GRAPH_HXX
#define ANDRES_GRAPH_HYPER_GRID_GRAPH_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <cstdlib>
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
    typedef std::vector<std::array<int, D>> OffsetVector;      //definition of OffsetVector
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
        size_type offsetIndex_;
    };

    /// AdjacencyIterator
    // \cond SUPPRESS_DOXYGEN
    class AdjacencyIterator
        : public std::iterator <
        std::random_access_iterator_tag,
        const AdjacencyType
        > {
    public:
        typedef HyperGridGraph<DIMENSION, Visitor> GraphType;
        typedef std::iterator <
            std::random_access_iterator_tag,
            const AdjacencyType
        > Base;
        typedef typename Base::difference_type difference_type;
        typedef typename Base::pointer pointer;
        typedef typename Base::reference reference;

        //constructors of AdjacencyIterator class
        AdjacencyIterator();
        AdjacencyIterator(const GraphType&);
        AdjacencyIterator(const GraphType&, const size_type);
        AdjacencyIterator(const GraphType&, const size_type, const size_type);

        // increment and decrement
        AdjacencyIterator& operator+=(const difference_type);
        AdjacencyIterator& operator-=(const difference_type);
        AdjacencyIterator& operator++(); // prefix
        AdjacencyIterator& operator--(); // prefix
        AdjacencyIterator operator++(int); // postfix
        AdjacencyIterator operator--(int); // postfix
        AdjacencyIterator operator+(const difference_type) const;
        AdjacencyIterator operator-(const difference_type) const;
        difference_type operator-(const AdjacencyIterator&) const;

        // comparison
        bool operator==(const AdjacencyIterator&) const;
        bool operator!=(const AdjacencyIterator&) const;
        bool operator<(const AdjacencyIterator&) const;
        bool operator<=(const AdjacencyIterator&) const;
        bool operator>(const AdjacencyIterator&) const;
        bool operator>=(const AdjacencyIterator&) const;
        
        // access
        reference operator*();
        pointer operator->();
        reference operator[](const difference_type);
        

    protected:
        
        const GraphType* graph_;
        //OffsetVector offsets_
        size_type vertex_;
        size_type adjacencyIndex_;
        AdjacencyType adjacency_;
    };

    //Defining the constructors of class HyperGridGraph
    HyperGridGraph(const Visitor & = Visitor());
    HyperGridGraph(const VertexCoordinate&, const OffsetVector&, const Visitor & = Visitor());
   // HyperGridGraph(const std::initializer_list<std::size_t>, const Visitor & = Visitor());//not implemented yet
    void assign(const Visitor & = Visitor());
    void assign(const VertexCoordinate&, const OffsetVector&, const Visitor & = Visitor());

    //iterator access
    AdjacencyIterator adjacenciesFromVertexBegin(const size_type) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const size_type) const;
    AdjacencyIterator adjacenciesToVertexBegin(const size_type) const;
    AdjacencyIterator adjacenciesToVertexEnd(const size_type) const;

    //access
    size_type numberOfVertices() const;
    size_type numberOfEdges() const;
    size_type numberOfEdgesFromVertex(const size_type) const;
    size_type numberOfEdgesToVertex(const size_type) const;
    size_type vertexFromVertex(const size_type, const size_type) const;
    size_type vertexToVertex(const size_type, const size_type) const;
    AdjacencyType adjacencyFromVertex(const size_type, const size_type) const;
    AdjacencyType adjacencyToVertex(const size_type, const size_type) const;
    std::pair<bool, size_type> findEdge(const size_type, const size_type) const;
    size_type insertEdge(const size_type, const size_type) const;

    size_type shape(const size_type) const;

    size_type vertex(const VertexCoordinate&) const;
    void vertex(size_type, VertexCoordinate&) const;
    size_type edge(const EdgeCoordinate&) const;
    void edge(size_type, EdgeCoordinate&) const;
    

private:

    size_type vertexFromVertex(const VertexCoordinate&,  const size_type, size_type&, bool&) const;
    void adjacencyFromVertex(const VertexCoordinate&,  const size_type, size_type&, size_type&) const;
    
    
    // Member variables
    VertexCoordinate shape_;
    OffsetVector offsets_;
    std::array<size_type, DIMENSION> edgeIndexOffsets_;
    std::array<size_type, DIMENSION> vertexIndexOffsets_;
    std::array<VertexCoordinate, DIMENSION> edgeShapes_; //in a HGG the value of 0 in a edgeShapes_ denotes that the edge is not possible for any rows
    size_type numberOfVertices_;
    Visitor visitor_;

};



/// Construct an empty hyper grid graph.
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
    : HyperGridGraph(VertexCoordinate({}), OffsetVector({}), visitor) // Chain-call Constructor
{}


/// Construct a hyper grid graph with a specified shape and offset.
///
/// \param shape the shape of the Hyper Grid Graph as an std::array.
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


// Clear a hyper grid graph. Assign function
///
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
///
template<unsigned char D, class VISITOR>
inline void
    HyperGridGraph<D, VISITOR>::assign(
        const Visitor& visitor
    ) {
    OffsetVector offset;
    VertexCoordinate shape;
    std::fill(shape.begin(), shape.end(), 0);
    offset.clear();
    assign(shape,offset, visitor);
}

/// Clear a hyper grid graph and assign a new shape.
///
/// \param shape the shape of the grid graph
/// \param 
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
 
template<unsigned char D, class VISITOR> 
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
    for (size_type i = 0; i < DIMENSION; ++i) {
        vertexIndexOffsets_[i] = cumprod;
        cumprod *= shape_[i];
    }
    numberOfVertices_ = cumprod; // Calculate the total number of vertices

    // Set edge offsets.
    size_type edgeIndexOffset = 0;
    for (size_type i = 0; i < offsets_.size(); ++i) {
        VertexCoordinate& edgeShape = edgeShapes_[i]; // Get the current edge shape array
        edgeShape = shape_; 
        
       
        for (size_type j = 0; j < DIMENSION; ++j) {
            int difference = edgeShape[j];
            difference = difference - offsets_[i][j];
            edgeShape[j] = abs(difference);
                   
        }
        // calculating cumulative product
        size_type cumprod = edgeShape[0];
        for (size_type k = 1; k < DIMENSION; ++k) {
            cumprod *= edgeShape[k];
        }
        edgeIndexOffsets_[i] = (edgeIndexOffset += cumprod);
    }       
   
}


//Calculation number of edges from vertex //completed
template<unsigned char D, class VISITOR> 
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::numberOfEdgesFromVertex(
    const size_type vertex
)   const {
    assert(vertex < numberOfVertices());
    VertexCoordinate edgeCoordinate;
    this->vertex(vertex, edgeCoordinate);
    size_type numEdgesFromVertex = 0;

    for (const auto& offset : offsets_)
    {
        VertexCoordinate neighborCoordinate = edgeCoordinate;
        VertexCoordinate negNeighborCoordinate = edgeCoordinate;
        bool isValidNeighbor1;
        bool isValidNeighbor2;

        // Add the elements of the current offset array to the edge coordinate  
        for (unsigned char i = 0; i < D; ++i)
        {
            neighborCoordinate[i] += offset[i];     
            // Check if the neighbor coordinate is within the shape of the graph
            isValidNeighbor1 = true;
       

            if (neighborCoordinate[i] < 0 || neighborCoordinate[i] >= shape_[i]) {
                isValidNeighbor1 = false;
                break;
            }
        }
        if (isValidNeighbor1) {
            ++numEdgesFromVertex;
        }

        for (unsigned char i = 0; i < D; ++i)
        {          
            negNeighborCoordinate[i] -= offset[i];

            // Check if the neighbor coordinate is within the shape of the graph
            isValidNeighbor2 = true;
            
            if (negNeighborCoordinate[i] < 0 || negNeighborCoordinate[i] >= shape_[i]) {
                isValidNeighbor2 = false;
                break;
            } 
        }
        if (isValidNeighbor2) {
            ++numEdgesFromVertex;
        }        
    }
    return numEdgesFromVertex;
}

/// Get the number of  edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::numberOfEdgesToVertex(
    const size_type vertex
) const {
    return numberOfEdgesFromVertex(vertex);
}

/// Get the integer index of a vertex from which a given vertex is reachable
/// via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and
/// numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::vertexFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    VertexCoordinate vertexCoordinate;
    this->vertex(vertex, vertexCoordinate);
    size_type direction;
    bool isSmaller;
    return vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
}


/// Get the integer index of a vertex to which a given vertex has a single
/// edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and
/// numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::vertexToVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return vertexFromVertex(vertex, j);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyType
HyperGridGraph<D, VISITOR>::adjacencyFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    size_type direction;
    size_type adjacentEdgeIndex;
    size_type adjacentVertexIndex;
    {
        VertexCoordinate vertexCoordinate;
        this->vertex(vertex, vertexCoordinate);

        adjacencyFromVertex(vertexCoordinate, j, adjacentVertexIndex, adjacentEdgeIndex);
    }
    return AdjacencyType(adjacentVertexIndex, adjacentEdgeIndex);
}


/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyType
HyperGridGraph<D, VISITOR>::adjacencyToVertex(
    const size_type vertex,
    const size_type j
) const {
    return adjacencyFromVertex(vertex, j);
}


/// Search for an edge (in constant time).
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \retval pair an \c std::pair. If an edge from \b vertex0 to \b vertex1
/// exists, \c pair.first is \c true
/// and \c pair.second is the index of such an edge. If no edge from \b
/// vertex0 to \b vertex1 exists, \c pair.first is \c false
/// and \c pair.second is undefined.
template<unsigned char D, class VISITOR>
inline std::pair<bool, typename HyperGridGraph<D, VISITOR>::size_type>
HyperGridGraph<D, VISITOR>::findEdge(
    const size_type vertex0,
    const size_type vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());
    VertexCoordinate vertexCoordinate0;
    vertex(vertex0, vertexCoordinate0);
    VertexCoordinate vertexCoordinate1;
    vertex(vertex1, vertexCoordinate1);
   
    for (size_type i = 0; i< offsets_.size(); i++) {
        bool flag1 = true;
        bool flag2 = true;
        for (size_t j = 0; j < DIMENSION; j++) {
            if (vertexCoordinate0[j] + offsets_[i][j] != vertexCoordinate1[j]) {
                flag1 = false;                
            }
            if (vertexCoordinate0[j] - offsets_[i][j] != vertexCoordinate1[j]) {
                flag2 = false;              
            }
        }
        if (flag1) {
            
            const EdgeCoordinate edgeCoordinate(vertexCoordinate0, i, false);
            const size_type edgeIndex = edge(edgeCoordinate);
            return  std::make_pair(true, edgeIndex);
        }
        if (flag2) {
       
            const EdgeCoordinate edgeCoordinate(vertexCoordinate1, i, false);
            const size_type edgeIndex = edge(edgeCoordinate);
            return  std::make_pair(true, edgeIndex);
        }        
    }
    return std::make_pair(false, 0);
}

/// Returns the edge between the specified vertices.
///
/// \param vertexIndex0 Integer index of the first vertex in the edge.
/// \param vertexIndex1 Integer index of the second vertex in the edge.
/// \return Integer index of the newly inserted edge.
/// \throw runtime_error If the edge does not exist.
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::insertEdge(
    const size_type vertexIndex0,
    const size_type vertexIndex1
) const {
    assert(vertexIndex0 < numberOfVertices());
    assert(vertexIndex1 < numberOfVertices());

    std::pair<bool, std::size_t> p = findEdge(vertexIndex0, vertexIndex1);
    if (p.first == true) {
        return p.second;
    }
    else {
        throw std::runtime_error("Edge not found.");
    }
}
/// Give Vertex index and retrieve the vertex coordinate. //completed
/// \param[in] vertexIndex the integer index of the requested vertex
/// \param[out] vertexCoordinate The coordinates of the vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<unsigned char D, class VISITOR> 
void
HyperGridGraph<D, VISITOR>::vertex(
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

/// Give Vertex Coordinate and retrieve the vertex index. //completed
/// \param vertexCoordinate the coordinate of the requested vertex
/// \return The integer index of the specified vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<unsigned char D, class VISITOR>
typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::vertex(
    const VertexCoordinate& vertexCoordinate
) const {
    size_type index = vertexCoordinate[DIMENSION - 1];
    for (size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * shape_[i - 1] + vertexCoordinate[i - 1];
    }
    return index;
}


// Get an iterator to the beginning of the sequence of adjacencies that
/// originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::adjacenciesFromVertexBegin(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies that originate
/// from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::adjacenciesFromVertexEnd(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies incident
/// to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::adjacenciesToVertexBegin(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies incident to
/// a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::adjacenciesToVertexEnd(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}



/// Get the number of vertices.

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::numberOfVertices(
) const {
return numberOfVertices_;
}

/// Get the number of edges.
///
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::numberOfEdges() const {
    return edgeIndexOffsets_[DIMENSION - 1];
}

/// Get the size of a specific dimension of the grid graph.
/// \param dimension the index of the dimension index to retrieve.
/// \return the size of the specified dimension.
/// \see mapIndexToCoordinate mapVertexCoordinateToIndex
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::shape(
    const size_type dimension
) const {
    return shape_[dimension];
}

/// Retrieve the specified edge of the graph. Given edge coordinate, return index
/// \param edgeCoordinate the coordinates of the minimum endpoint (\e pivot)
/// and the direction of the requested edge.
/// \return The integer index of the specified edge.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa hasEdge()
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::edge(
    const EdgeCoordinate& edgeCoordinate
) const {
    assert(edgeCoordinate.offsetIndex_ < offsets_.size());

    const size_type& offsetIndex_ = edgeCoordinate.offsetIndex_;
    const VertexCoordinate& pivotCoordinate = edgeCoordinate.pivotCoordinate_;
    const VertexCoordinate& edgeShape = edgeShapes_[offsetIndex_];

    size_type index = pivotCoordinate[DIMENSION - 1];
    for (size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * edgeShape[i - 1] + pivotCoordinate[i - 1];
    }
    if (offsetIndex_ > 0) {
        const size_type& indexOffset = edgeIndexOffsets_[offsetIndex_ - 1];
        index += indexOffset;
    }
    return index;
}



/// Retrieve the specified edge of the graph. Given edge index, send back edgeCoordinate(Pivot coordinate and offset index).
/// \param[in] edgeIndex the integer index of the requested edge.
/// \param[out] edgeCoordinate a HyperGridGraph::EdgeCoordinate instance. \c
/// edgeCoordinate.pivot is the integer index
///  of the minimum of the two edge endpoints; \p edgeCoordinate.direction
///  is the dimension along which
///  the edge is drawn, with an assumed positive isSmaller.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa HyperGridGraph::bool
template<unsigned char D, class VISITOR>
inline void
HyperGridGraph<D, VISITOR>::edge(
        size_type edgeIndex,
        EdgeCoordinate & edgeCoordinate
    ) const {
    // WARNING: If the assertion does not hold, the code will scan until an unlikely condition.
    assert(edgeIndex < numberOfEdges());

    size_type& direction = edgeCoordinate.offsetIndex_;
    // Find the direction as the last edge offset:
    for (direction = 0; edgeIndex >= edgeIndexOffsets_[direction]; ++direction);
    if (direction > 0) { // Not needed, but saves one memory lookup.
        const size_type& offset = edgeIndexOffsets_[direction - 1];
        edgeIndex -= offset;
    }
    const VertexCoordinate& edgeShape = edgeShapes_[direction];
    // mapEdgeIndexToCoordinate
    {
        VertexCoordinate& pivotCoordinate = edgeCoordinate.pivotCoordinate_;
        size_type i;
        for (i = 0; i < DIMENSION - 1; ++i) {
            pivotCoordinate[i] = edgeIndex % edgeShape[i];
            edgeIndex = edgeIndex / edgeShape[i];
        }
        pivotCoordinate[i] = edgeIndex;
    }
}

/// \internal
/// \brief Retrieve the direction and isSmaller of the <b>j</b>-th edge of
/// a vertex
/// \param vertexCoordinate the integer index of the vertex from which the
/// edges are originating
/// \param j the index within the set of adacent edges of the specified vertex
/// \param[out] direction the direction of the specified edge
/// \param[out] isSmaller the isSmaller of the specified edge
/// \retval found If the edge is found, a value of \c true is returned.
/// \remark This function attempts a fast imlementation that tries not
/// to waste any comparison operations and is meant for use as a building
/// block of the class.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa directionOfEdgeFromVertex, edgeFromVertex
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::size_type
HyperGridGraph<D, VISITOR>::vertexFromVertex(
    const VertexCoordinate& vertexCoordinate,
    const size_type j,
    size_type& direction,
    bool& isSmaller
) const {
    assert(vertex(vertexCoordinate) < numberOfVertices());
    VertexCoordinate modifieableVertexCoordinate = vertexCoordinate;
    size_type cur_j = 0;
   
    for (size_type i = 0; i < offsets_.size(); i++) //incoming edges
    {
        bool flag = true;
        for (size_type k = DIMENSION; k > 0; k--)
        {
            int z = (int)vertexCoordinate[k - 1] - (int)offsets_[i][k - 1];
            if (z < 0) //check if edge doesn't exist in this dimension
            {
                flag = false; 
                break;
            }
        }
        if (flag == true && cur_j == j)
        {
            for (size_type k = 0; k < DIMENSION; k++)
            {
                modifieableVertexCoordinate[k] -= offsets_[i][k];  //modify the vertex
            }
            direction = i; //along i-th offset
            isSmaller = false;
            const size_type vertexIndex = vertex(modifieableVertexCoordinate);
            return vertexIndex;
        }
        if (flag)
        {
            ++cur_j;
        }        
    }

    for (size_type i = 0; i < offsets_.size(); i++) // outgoing edges
    {
        bool flag = true;
        for (size_type k = 0; k < DIMENSION; k++)
        {
            int z = (int)vertexCoordinate[k] + (int)offsets_[i][k];
            if (z >= shape_[k])
            {
                flag = false;
                break;
            }
        }
        if (flag == true  && cur_j == j)
        {
            for (size_type k = 0; k < DIMENSION; k++)
            {
                modifieableVertexCoordinate[k] += offsets_[i][k]; //modify the vertex
            }
            direction = i; //along i-th offset
            isSmaller = true;
            const size_type vertexIndex = vertex(modifieableVertexCoordinate);
            return vertexIndex;
        }
        if (flag)
        {
            ++cur_j;
        }               
    }
    throw std::out_of_range("vertex neighbor index out of range.");
}


template<unsigned char D, class VISITOR>
inline void
HyperGridGraph<D, VISITOR>::adjacencyFromVertex(
    const VertexCoordinate& vertexCoordinate,
    const size_type j,
    size_type& adjacentVertexIndex,
    size_type& adjacentEdgeIndex
) const {
    assert(vertex(vertexCoordinate) < numberOfVertices());
    size_type direction;
    bool isSmaller;
    adjacentVertexIndex = vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
    if (isSmaller) //outgoing edge
    {
        
        adjacentEdgeIndex = edge(EdgeCoordinate(vertexCoordinate, direction));
    }
    else { //incoming edge
        VertexCoordinate pivotVertexCoordinate;
        this->vertex(adjacentVertexIndex, pivotVertexCoordinate);
        adjacentEdgeIndex = edge(EdgeCoordinate(pivotVertexCoordinate, direction));
    }
}



/// Initialize an edge coordinate.
/// \param pivotCoordinate coordinate of the reference vertex.
/// \param offsetIndex_ the index of the offset along which the edge is drawn.
/// \param isSmaller relative position of specified endpoint.\n
/// A value of \c true results in an object corresponding to the edge
/// of which the smallest endpoint has a HyperGridGraph::VertexCoordinate
/// equal to \b pivotCoordinate.\n
/// Similarly, a value of \c false results in an object corresponding
/// to the edge of which the largest endpoint has a
/// GridGraph::VertexCoordinate equal to \b pivotCoordinate.
template<unsigned char D, class VISITOR>
inline HyperGridGraph<D, VISITOR>::EdgeCoordinate::EdgeCoordinate(
    const VertexCoordinate& pivotCoordinate,
    const size_type offsetIndex_,
    const bool isSmaller
)
    : pivotCoordinate_(pivotCoordinate),
    offsetIndex_(offsetIndex_) {
    if (isSmaller) {
        assert(pivotCoordinate_[offsetIndex_] > 0);
        --pivotCoordinate_[offsetIndex_];
    }
}

/// Default non-initializing constructor
/// \warning this constructor will \b NOT set the pivotCoordinate_ membr
/// variable to zero.
template<unsigned char D, class VISITOR>
inline HyperGridGraph<D, VISITOR>::EdgeCoordinate::EdgeCoordinate() {
}





//Adjacency Iterator implementation

template<unsigned char D, class VISITOR>
inline
HyperGridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator()
    : vertex_(0),
    adjacencyIndex_(0),
    adjacency_()
{}
template<unsigned char D, class VISITOR>
inline
HyperGridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph
)
    : graph_(&graph),
    vertex_(0),
    adjacencyIndex_(0),
    adjacency_()
{}

template<unsigned char D, class VISITOR>
inline
HyperGridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const size_type vertex
)
    : graph_(&graph),
    vertex_(vertex),
    adjacencyIndex_(0),
    adjacency_() {
    assert(vertex < graph.numberOfVertices());
}

template<unsigned char D, class VISITOR>
inline
HyperGridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
    : graph_(&graph),
    vertex_(vertex),
    adjacencyIndex_(adjacencyIndex),
    adjacency_() {
    assert(vertex < graph.numberOfVertices());
    assert(adjacencyIndex <= graph.numberOfVertices());
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator&
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator+=(
    const difference_type d
    ) {
    adjacencyIndex_ += d;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator&
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator-=(
    const difference_type d
    ) {
    adjacencyIndex_ -= d;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator&
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator++() {
    ++adjacencyIndex_;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator&
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator--() {
    --adjacencyIndex_;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator++(int) {
    AdjacencyIterator copy = *this;
    ++adjacencyIndex_;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator--(int) {
    AdjacencyIterator copy = *this;
    --adjacencyIndex_;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator+(
    const difference_type d
    ) const {
    AdjacencyIterator copy = *this;
    copy += d;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator-(
    const difference_type d
    ) const {
    AdjacencyIterator copy = *this;
    copy -= d;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator::difference_type
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator-(
    const AdjacencyIterator& adjacencyIterator
    ) const {
    return adjacencyIndex_ - adjacencyIterator.adjacencyIndex_;
}

template<unsigned char D, class VISITOR>
inline bool
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator==(
    const AdjacencyIterator& other
    ) const {
    return adjacencyIndex_ == other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator!=(
    const AdjacencyIterator& other
    ) const {
    return adjacencyIndex_ != other.adjacencyIndex_
        || vertex_ != other.vertex_
        || graph_ != other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator<(
    const AdjacencyIterator& other
    ) const {
    return adjacencyIndex_ < other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator<=(
    const AdjacencyIterator& other
    ) const {
    return adjacencyIndex_ <= other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator>(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ > other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator>=(
    const AdjacencyIterator& other
    ) const {
    return adjacencyIndex_ >= other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}


// Since the HyperGridGraph has no backing storage for each of the Adjacency
// objects handled in the class, the AdjacencyIterator is limited to only
// one adjacency object per iterator instance.
// This should be sufficient for most normal usage scenarios, and is supported
// by the very lightweight copying of the AdjacencyType object itself.
// Note, however, that if you store a reference to the pointed object,
// then advance the iterator and subsequently dereference it again,
// the new reference will refer to the very same same adjacency object,
// unique for the AdjacencyIterator instance.
// This operation will therefore silently update all previous references to
// the adjacency object of the iterator.
template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator::reference
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator*() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return adjacency_;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator::pointer
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator->() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return &adjacency_;
}

template<unsigned char D, class VISITOR>
inline typename HyperGridGraph<D, VISITOR>::AdjacencyIterator::reference
HyperGridGraph<D, VISITOR>::AdjacencyIterator::operator[](
    const difference_type j) {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_ + j);
    return adjacency_;
}


} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_HYPER_GRID_GRAPH_HXX