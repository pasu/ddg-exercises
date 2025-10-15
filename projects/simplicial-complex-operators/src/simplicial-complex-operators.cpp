// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    std::vector<Eigen::Triplet<size_t>> triplets;
    for (Edge e : mesh->edges()) {
        // For each vertex, iterate over its incident edges
        for (Vertex v : e.adjacentVertices()) {
            // Add an entry to the adjacency matrix
            Eigen::Triplet<size_t> ele(e.getIndex(), v.getIndex(), 1);
            triplets.push_back(ele);
        }
    }

    SparseMatrix<size_t> A_Vextex_Edge(mesh->nEdges(), mesh->nVertices());
    A_Vextex_Edge.setFromTriplets(triplets.begin(), triplets.end());
    return A_Vextex_Edge;

    // return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    std::vector<Eigen::Triplet<size_t>> triplets;
    for (Face f : mesh->faces()) {
        // For each face, iterate over its incident edges
        for (Edge e : f.adjacentEdges()) {
            // Add an entry to the adjacency matrix
            Eigen::Triplet<size_t> ele(f.getIndex(), e.getIndex(), 1);
            triplets.push_back(ele);
        }
    }

    SparseMatrix<size_t> A_Edge_Face(mesh->nFaces(), mesh->nEdges());
    A_Edge_Face.setFromTriplets(triplets.begin(), triplets.end());
    return A_Edge_Face;
    // return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vertexVector = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t i : subset.vertices) {
        vertexVector(i) = 1;
    }
    return vertexVector;
    // return Vector<size_t>::Zero(1);
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> edgeVector = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t i : subset.edges) {
               edgeVector(i) = 1;
    }
    return edgeVector;
    // return Vector<size_t>::Zero(1);
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> faceVector = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t i : subset.faces) {
               faceVector(i) = 1;
    }
    return faceVector;
    // return Vector<size_t>::Zero(1);
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset) + A0 * vertexVector;
    Vector<size_t> faceVector = buildFaceVector(subset) + A1 * edgeVector;

    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    for (size_t i = 0; i < vertexVector.size(); i++) {
        if (vertexVector(i) > 0) {
            vertices.insert(i);
        }
    }

    for (size_t i = 0; i < edgeVector.size(); i++) {
        if (edgeVector(i) > 0) {
            edges.insert(i);
        }
    }

    for (size_t i = 0; i < faceVector.size(); i++) {
        if (faceVector(i) > 0) {
            faces.insert(i);
        }
    }

    return MeshSubset(vertices, edges, faces);
    // return subset; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset) + A1.transpose() * faceVector;
    Vector<size_t> vertexVector = buildVertexVector(subset) + A0.transpose() * edgeVector;

    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    for (size_t i = 0; i < vertexVector.size(); i++) {
        if (vertexVector(i) > 0) {
            vertices.insert(i);
        }
    }

    for (size_t i = 0; i < edgeVector.size(); i++) {
        if (edgeVector(i) > 0) {
            edges.insert(i);
        }
    }

    for (size_t i = 0; i < faceVector.size(); i++) {
        if (faceVector(i) > 0) {
            faces.insert(i);
        }
    }

    return MeshSubset(vertices, edges, faces);
    // return subset; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset closureSubset = closure(star(subset));
    MeshSubset starSubset = star(closure(subset));

    closureSubset.deleteSubset(starSubset);

    return closureSubset;

    // return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return subset.equals(closure(subset));
    // return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO

    // Early exit if the subset is not a valid complex
    if (!isComplex(subset)) {
        return -1;
    }

    // Build vectors for simplices
    auto faceA = buildFaceVector(subset);
    auto edgeA = buildEdgeVector(subset);
    auto vertexA = buildVertexVector(subset);

    Vector<size_t> edgeB = (faceA.transpose() * A1).cwiseMin(1);
    Vector<size_t> vertexB = (edgeA.transpose() * A0).cwiseMin(1);

    /* Compute expected edges and vertices from faces and edges
    auto edgeB = (faceA.transpose() * A1).cwiseMin(1);
    auto vertexB = (edgeA.transpose() * A0).cwiseMin(1);*/


    if (faceA.sum() > 0) {
        return edgeA == edgeB && vertexA == vertexB ? 2 : -1;
    }

    if (edgeA.sum() > 0) {
        return vertexA == vertexB ? 1 : -1;
    }

    if (vertexA.sum() > 0) {
        return 0;
    }

    return -1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    // TODO

    // Compute edge and vertex vectors
    Vector<size_t> edgeVector = A1.transpose() * buildFaceVector(subset);
    Vector<size_t> vertexVector = A0.transpose() * buildEdgeVector(subset);

    // Initialize the boundary subset
    MeshSubset boundarySubset;

    // Add edges to the boundary subset where edgeVector[i] == 1
    for (size_t edgeIndex = 0; edgeIndex < mesh->nEdges(); ++edgeIndex) {
        if (edgeVector[edgeIndex] == 1) {
            boundarySubset.addEdge(edgeIndex);
        }
    }

    // Add vertices to the boundary subset where vertexVector[i] == 1
    for (size_t vertexIndex = 0; vertexIndex < mesh->nVertices(); ++vertexIndex) {
        if (vertexVector[vertexIndex] == 1) {
            boundarySubset.addVertex(vertexIndex);
        }
    }

    // Return the closure of the boundary subset
    return closure(boundarySubset);

    // return subset; // placeholder
}