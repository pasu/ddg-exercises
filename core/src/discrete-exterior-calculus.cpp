// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mesh.nVertices());

    for (Vertex v : mesh.vertices()) {
		double area = vertexDualArea(v);
		tripletList.push_back(Eigen::Triplet<double>(v.getIndex(), v.getIndex(), area));
	}

    SparseMatrix<double> H0(mesh.nVertices(), mesh.nVertices());
    H0.setFromTriplets(tripletList.begin(), tripletList.end());
    return H0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mesh.nEdges());
    for (Edge e : mesh.edges()) {
        double area = (cotan(e.halfedge()) + cotan(e.halfedge().twin())) / 2.0;
        tripletList.push_back(Eigen::Triplet<double>(e.getIndex(), e.getIndex(), area));
    }

    SparseMatrix<double> H1(mesh.nEdges(), mesh.nEdges());
    H1.setFromTriplets(tripletList.begin(), tripletList.end());
    return H1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mesh.nFaces());
    for (Face f : mesh.faces()) {
		double area = 1.0 / faceArea(f);
		tripletList.push_back(Eigen::Triplet<double>(f.getIndex(), f.getIndex(), area));
	}
    SparseMatrix<double> H2(mesh.nFaces(), mesh.nFaces());
    H2.setFromTriplets(tripletList.begin(), tripletList.end());
    return H2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mesh.nEdges() * 2);
    for (Edge e : mesh.edges()) {
		Halfedge he = e.halfedge();
		Vertex vStart = he.vertex();
		Vertex vEnd = he.next().vertex();
		tripletList.push_back(Eigen::Triplet<double>(e.getIndex(), vStart.getIndex(), -1.0));
		tripletList.push_back(Eigen::Triplet<double>(e.getIndex(), vEnd.getIndex(), 1.0));
	}

    SparseMatrix<double> D0(mesh.nEdges(), mesh.nVertices());
    D0.setFromTriplets(tripletList.begin(), tripletList.end());
    return D0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mesh.nFaces() * 3);
    for (Face f : mesh.faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            double orientation = he.getIndex() < he.twin().getIndex() ? 1 : -1;
            tripletList.push_back(Eigen::Triplet<double>(f.getIndex(), he.edge().getIndex(), orientation));
        } 
    }
    SparseMatrix<double> D1(mesh.nFaces(), mesh.nEdges());
    D1.setFromTriplets(tripletList.begin(), tripletList.end());
    return D1;
}

} // namespace surface
} // namespace geometrycentral