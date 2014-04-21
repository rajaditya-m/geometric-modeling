// ------------------------------------------------------------
// Mesh.h
// ------------------------------------------------------------
//
// ------------------------------------------------------------

#ifndef MESHCLASS
#define	MESHCLASS

#include <stdio.h>
#include <vector>
#include <set>


// Change to double for more precision
typedef		float	datatype;


// ------------------------------------------------------------
// GeomVert: this class holds the geometric coordinates of a vertex
// ------------------------------------------------------------
class GeomVert {
public:
	GeomVert(datatype x, datatype y, datatype z) { mCo[0] = x; mCo[1] = y; mCo[2] = z; }

	datatype      GetCo(int axis)                { return mCo[axis]; }

	bool operator == (GeomVert &A)               { 
		return ( (mCo[0] == A.GetCo(0)) && (mCo[1] == A.GetCo(1)) && (mCo[2] == A.GetCo(2)) );		
	}
 
private:
	datatype	mCo[3];
};
// ------------------------------------------------------------





// ------------------------------------------------------------
// TopoVert: this class holds all of a vertex's topological (connectivity) information
// ------------------------------------------------------------
class TopoVert {
public:
	TopoVert()                      { };
	~TopoVert()                     { mIncVerts.clear(); mIncEdges.clear(); mIncFacets.clear(); }
	void AddIncVert (int vert_ind ) { mIncVerts.insert( vert_ind ); }
	//void AddIncVert (int vert_ind ) { mIncVerts.push_back( vert_ind ); }
	void AddIncEdge (int edge_ind ) { mIncEdges.push_back( edge_ind ); }
	void AddIncFacet(int facet_ind) { mIncFacets.push_back( facet_ind ); }
	
	int GetNumberIncVertices()		{ return mIncVerts.size(); }
	int GetIncVertex(int vert_ind)  { 
		std::set<int>::iterator sit = mIncVerts.begin(); for (int i = 0; i < vert_ind; i++) sit++;		
		return *sit;
	}
	int GetNumberIncEdges()			{ return mIncEdges.size(); }
	int GetIncEdge(int edge_ind)    { return mIncEdges[edge_ind]; }
	int GetNumberIncFacets()	    { return mIncFacets.size(); }
	int GetIncFacet(int facet_ind)  { return mIncFacets[facet_ind]; }

private:
	std::set<int>    mIncVerts;
	std::vector<int> mIncEdges;
	std::vector<int> mIncFacets;
	
};
// ------------------------------------------------------------





// ------------------------------------------------------------
// TopoEdge
// ------------------------------------------------------------
class TopoEdge {
public:
	TopoEdge()                      { v1 = v2 = -1; }
	~TopoEdge()                     { mIncFacets.clear(); }

	bool operator == (TopoEdge &A)  {
		return (  ((v1 == A.GetVertex(0)) && (v2 == A.GetVertex(1))) || ((v2 == A.GetVertex(0)) && (v1 == A.GetVertex(1))) );
	}

	int  GetVertex(int ind)         { if (ind == 0) return v1;  return v2; }
	void SetVertex(int ind, int v)  { if (ind == 0) { v1 = v; } else { v2 = v; } }
	
	void AddIncFacet(int facet_ind) { mIncFacets.push_back(facet_ind); }
	int  GetNumberIncFacets()       { return mIncFacets.size(); }
	int  GetIncFacet(int facet_ind) { return mIncFacets[facet_ind]; }

private:
	int v1, v2;
	std::vector<int> mIncFacets;
};
// ------------------------------------------------------------





// ------------------------------------------------------------
// TopoFacet:  this class holds a facet's topological connectivity) information
//             Facets are represented as a list of vertex indices
// ------------------------------------------------------------
class TopoFacet {
public:
	TopoFacet()                     { };
	~TopoFacet()                    { mIncVerts.clear(); mIncEdges.clear();  mIncFacets.clear(); }
	void AddIncVertex(int v_ind)    { mIncVerts.push_back( v_ind ); }
	void AddIncEdge(int e_ind)      { mIncEdges.push_back( e_ind ); }
	void AddIncFacet(int f_ind)     { mIncFacets.insert( f_ind ); }
	int  GetNumberVertices()        { return mIncVerts.size(); }
	int  GetVertexInd(int vert_ind) { return mIncVerts[vert_ind]; }
	int  GetNumberEdges()		    { return mIncEdges.size(); }
	int  GetIncEdge(int edge_ind)   { return mIncEdges[edge_ind]; }
	int  GetNumberFacets()		    { return mIncFacets.size(); }
	int  GetIncFacet(int facet_ind) { 
		std::set<int>::iterator sit = mIncFacets.begin(); for (int i = 0; i < facet_ind; i++) sit++;		
		return *sit;
	}



private:
	std::vector<int> mIncVerts;	
	std::vector<int> mIncEdges;
	std::set<int>    mIncFacets;
};
// ------------------------------------------------------------






// ------------------------------------------------------------
// Mesh:  This class uses all the preceding classes to represent a mesh with
//        adjacency.connectivity information
// ------------------------------------------------------------
class Mesh {
public:
	Mesh()  { };
	~Mesh() { Erase(); };

	void      AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2, datatype x3, datatype y3, datatype z3);
	void      AddFacet(std::vector<GeomVert> geomfacet);
	
	int		  GetNumberVertices()		  { return mGeomVerts.size(); }
	int		  GetNumberEdges()			  { return mTopoEdges.size(); }
	int       GetNumberFacets()           { return mTopoFacets.size(); }

	TopoVert  GetVertex(int vert_ind)     { return mTopoVerts[vert_ind]; }
	TopoEdge  GetEdge(int edge_ind)       { return mTopoEdges[edge_ind]; }
	TopoFacet GetFacet(int facet_ind)     { return mTopoFacets[facet_ind]; }
	
	GeomVert  GetGeomVertex(int vert_ind) { return mGeomVerts[vert_ind]; }



private:
	int       FindGeomVertex(GeomVert v);
	int		  FindTopoEdge(TopoEdge e);
	void      Erase();

	std::vector<GeomVert>  mGeomVerts;
	std::vector<TopoVert>  mTopoVerts;
	std::vector<TopoEdge>  mTopoEdges;
	std::vector<TopoFacet> mTopoFacets;
};
// ------------------------------------------------------------

#endif