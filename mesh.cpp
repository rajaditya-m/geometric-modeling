// ------------------------------------------------------------
// Mesh.cpp:  Implementation of mesh class
// ------------------------------------------------------------

#include "mesh.h"





// ------------------------------------------------------------
// AddFacet:  Adds a triangle to the mesh.
//            This is one of 2 functions that can be used to build a mesh
// ------------------------------------------------------------
void Mesh::AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2, datatype x3, datatype y3, datatype z3) {
	
	std::vector<GeomVert> geomfacet;
	geomfacet.push_back( GeomVert(x1, y1, z1) );
	geomfacet.push_back( GeomVert(x2, y2, z2) );
	geomfacet.push_back( GeomVert(x3, y3, z3) );

	AddFacet( geomfacet );
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// AddFacet:  Adds a facet with arbitrary number of vertices to mesh
//            This is one of 2 functions that can be used to build a mesh
// ------------------------------------------------------------
void Mesh::AddFacet(std::vector<GeomVert> geomfacet) {
	int i;

	// --------------
	// Create topo facet (list of geom vertex indices)
	TopoFacet topofacet;
	// Look for facet vertices in mesh - if they don't already exist in mesh then add them
	for (i = 0; i < geomfacet.size(); i++) {
		int v_ind = FindGeomVertex( geomfacet[i] );
		if (v_ind == -1) {
			// New vertex:  add geomtric vertex
			v_ind = mGeomVerts.size();
			mGeomVerts.push_back( geomfacet[i] );

			// Add topo vertex
			TopoVert topovert;
			mTopoVerts.push_back( topovert );
		}

		// Add vertex indice to topo facet
		topofacet.AddIncVertex( v_ind );
	}

	// Add this new topo facet to mesh	
	int facet_ind = mTopoFacets.size();
	mTopoFacets.push_back( topofacet );


	// Add edges of facet to mesh, again checking if they already exist
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1  :  i - 1;
		
		// Create edge
		TopoEdge e;
		e.SetVertex(0, topofacet.GetVertexInd(prev) );
		e.SetVertex(1, topofacet.GetVertexInd(i) );

		// Check if exists
		int e_ind = FindTopoEdge( e );
		
		if (e_ind == -1) {
			// Didn't exist, add to mesh
			e_ind = mTopoEdges.size();
			mTopoVerts[ e.GetVertex(0) ].AddIncEdge( e_ind );
			mTopoVerts[ e.GetVertex(1) ].AddIncEdge( e_ind );
			mTopoEdges.push_back( e );			
		}

		// Point edge to this facet
		mTopoEdges[e_ind].AddIncFacet( facet_ind );

		// Point facet to this edge
		mTopoFacets[ facet_ind ].AddIncEdge( e_ind );
	}
	// --------------
		

	
	// Compute other connectivity
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		// Add vertex-facet topology
		mTopoVerts[  topofacet.GetVertexInd(i) ].AddIncFacet( facet_ind );

		// Add vertex-vertex (edge) topology
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1  :  i - 1;
		int next = (i == topofacet.GetNumberVertices() - 1) ? 0 : i + 1;

		mTopoVerts[ topofacet.GetVertexInd(i) ].AddIncVert( topofacet.GetVertexInd( prev ) );
		mTopoVerts[ topofacet.GetVertexInd(i) ].AddIncVert( topofacet.GetVertexInd( next ) );
	}
	
	// Facet-facet adjacency...
	for (i = 0; i < mTopoFacets[ facet_ind ].GetNumberEdges(); i++) {		
		TopoEdge edge = mTopoEdges[ mTopoFacets[ facet_ind ].GetIncEdge(i) ];
		for (int j = 0; j < edge.GetNumberIncFacets(); j++) {
			if (edge.GetIncFacet(j) != facet_ind) {
				mTopoFacets[ facet_ind ].AddIncFacet( edge.GetIncFacet(j) );
				mTopoFacets[ edge.GetIncFacet(j) ].AddIncFacet( facet_ind );
			}
		}
	}
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// Erase:  Releases all memory used by object
// ------------------------------------------------------------
void Mesh::Erase() {
	mGeomVerts.clear();
	mTopoVerts.clear();
	mTopoEdges.clear();
	mTopoFacets.clear();
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// FindGeomVertex:  Searches for a geometric vertex in the mesh,
//                  returning its indice if found, -1 otherwise
// ------------------------------------------------------------
int Mesh::FindGeomVertex(GeomVert v) {
	for (int i = 0; i < mGeomVerts.size(); i++) {
		if (mGeomVerts[i] == v) return i;
	}
	return -1;
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// FindTopoEdge:  Searches for an edge in the mesh, returing 
//                its indice if found, -1 otherwise
// ------------------------------------------------------------
int	Mesh::FindTopoEdge(TopoEdge e) {
	for (int i = 0; i < mTopoEdges.size(); i++) {
		if (mTopoEdges[i] == e) return i;
	}
	return -1;
}
// ------------------------------------------------------------

