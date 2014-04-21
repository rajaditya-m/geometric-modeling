#ifndef __SUBDIVALGO_H_
#define __SUBDIVALGO_H_

#include "mesh.h"

#define ALPHA 0.625f

glm::vec3 doo_sabine_vf(int vertex_id,int facet_id,Mesh* mesh)
{
	//Compute centroid of the face 
	glm::vec3 x(0.0,0.0,0.0);
	TopoFacet current_face = mesh->GetFacet(facet_id);
	int n_vertices = current_face.GetNumberVertices();
	for(int i=0;i<n_vertices;i++)
	{
		int cur_v_idx = current_face.GetVertexInd(i);
		GeomVert cur_v = mesh->GetGeomVertex(cur_v_idx);
		glm::vec3 p(cur_v.GetCo(0),cur_v.GetCo(1),cur_v.GetCo(2));
		x += p;
	}
	float avg_fac = 1.0f/(float)n_vertices ; 
	x = avg_fac*x;

	//Get the two face centroid that is ve1 and ve2
	int ve1_idx,ve2_idx;
	bool ve1_found = false;
	bool ve2_found = false;
	int n_edges = current_face.GetNumberEdges();
	for(int i=0;i<n_edges;i++)
	{
		int cur_e_idx = current_face.GetIncEdge(i);
		TopoEdge cur_e = mesh->GetEdge(cur_e_idx);
		if((cur_e.GetVertex(0)==vertex_id)||(cur_e.GetVertex(1)==vertex_id))
		{
			int n_surrounding_faces = cur_e.GetNumberIncFacets();
			for(int k=0;k<n_surrounding_faces;k++)
			{
				int surr_f_e = cur_e.GetIncFacet(k);
				if(surr_f_e==facet_id)
				{
					if(!ve1_found)
					{
						ve1_idx = cur_e_idx;
						ve1_found = true;
						break; //Inner loop
					}
					else
					{
						ve2_idx = cur_e_idx;
						ve2_found = true;
						break; //Inner loop
					}
				}
			}
		}
		if(ve2_found)
			break;
	} 
	TopoEdge e1 = mesh->GetEdge(ve1_idx);
	TopoEdge e2 = mesh->GetEdge(ve2_idx);
	GeomVert e1_1 = mesh->GetGeomVertex(e1.GetVertex(0));
	GeomVert e1_2 = mesh->GetGeomVertex(e1.GetVertex(1));
	GeomVert e2_1 = mesh->GetGeomVertex(e2.GetVertex(0));
	GeomVert e2_2 = mesh->GetGeomVertex(e2.GetVertex(1));
	glm::vec3 e11(e1_1.GetCo(0),e1_1.GetCo(1),e1_1.GetCo(2));
	glm::vec3 e12(e1_2.GetCo(0),e1_2.GetCo(1),e1_2.GetCo(2));
	glm::vec3 e21(e2_1.GetCo(0),e2_1.GetCo(1),e2_1.GetCo(2));
	glm::vec3 e22(e2_2.GetCo(0),e2_2.GetCo(1),e2_2.GetCo(2));
	glm::vec3 ve1 = 0.5f*(e11+e12);
	glm::vec3 ve2 = 0.5f*(e21+e22);

	GeomVert vg = mesh->GetGeomVertex(vertex_id);
	glm::vec3 v(vg.GetCo(0),vg.GetCo(1),vg.GetCo(2));

	glm::vec3 res = 0.25f*(x+v+ve1+ve2);

	return res;
}

void doo_sabine(Mesh* input,Mesh** output)
{
	//Create face-face
	int n_faces = input->GetNumberFacets();
	for(int i=0;i<n_faces;i++)
	{
		TopoFacet face = input->GetFacet(i);
		std::vector<GeomVert> temp_face_buffer;
		int n_vertices_in_face = face.GetNumberVertices();
		for(int j=0;j<n_vertices_in_face;j++)
		{
			int vert_id = face.GetVertexInd(j);
			glm::vec3 pnt = doo_sabine_vf(vert_id,i,input);
			temp_face_buffer.push_back(GeomVert(pnt.x,pnt.y,pnt.z));
		}
		(*output)->AddFacet(temp_face_buffer);
	}
	//Create edge-face
	int n_edges = input->GetNumberEdges();
	for(int i=0;i<n_edges;i++)
	{
		TopoEdge edge = input->GetEdge(i);
		std::vector<GeomVert> temp_face_buffer;
		int n_faces_surrounding = edge.GetNumberIncFacets();
		if(n_faces_surrounding==1)
			continue;
		else if(n_faces_surrounding==2)
		{
			int face_id_1 = edge.GetIncFacet(0);
			int face_id_2 = edge.GetIncFacet(1);
			int vert_id_1 = edge.GetVertex(0);
			int vert_id_2 = edge.GetVertex(1);
			glm::vec3 v1f1 = doo_sabine_vf(vert_id_1,face_id_1,input);
			temp_face_buffer.push_back(GeomVert(v1f1.x,v1f1.y,v1f1.z));
			glm::vec3 v1f2 = doo_sabine_vf(vert_id_1,face_id_2,input);
			temp_face_buffer.push_back(GeomVert(v1f2.x,v1f2.y,v1f2.z));
			glm::vec3 v2f2 = doo_sabine_vf(vert_id_2,face_id_2,input);
			temp_face_buffer.push_back(GeomVert(v2f2.x,v2f2.y,v2f2.z));
			glm::vec3 v2f1 = doo_sabine_vf(vert_id_2,face_id_1,input);
			temp_face_buffer.push_back(GeomVert(v2f1.x,v2f1.y,v2f1.z));
		}
		else
		{
			std::cout << "Something is wrong here.";
		}
		(*output)->AddFacet(temp_face_buffer);
	}
	//Create vertex-face
	int n_vertices = input->GetNumberVertices();
	for(int i=0;i<n_vertices;i++)
	{
		std::set<int> already_processed;
		TopoVert ver = input->GetVertex(i);
		std::vector<GeomVert> temp_face_buffer;
		int n_faces_surrounding = ver.GetNumberIncFacets();
		int seed =  ver.GetIncFacet(0);
		bool next_seed_found = false;
		// std::cout << "Face indices are " ;
		// for(int j=0;j<n_faces_surrounding;j++)
		// {
		// 	std::cout << ver.GetIncFacet(j) << " ";
		// }
		// std::cout << "\n";
		// std::cout << "Seed indices are " ;
		for(int j=0;j<n_faces_surrounding;j++)
		{
			already_processed.insert(seed);
			// std::cout <<  seed << " ";
			next_seed_found = false;
			// int face_id = ver.GetIncFacet(seed);
			glm::vec3 pnt = doo_sabine_vf(i,seed,input);
			temp_face_buffer.push_back(GeomVert(pnt.x,pnt.y,pnt.z));
			//Generate the next seed
			//Look at all the surrounding faces of this face
			TopoFacet cur_seed = input->GetFacet(seed);
			int n_sur_face_cur_face = cur_seed.GetNumberFacets();
			for(int k=0;k<n_sur_face_cur_face;k++)
			{
				int idx_facet = cur_seed.GetIncFacet(k);
				std::set<int>::iterator it;
				it = already_processed.find(idx_facet);
				if(it!=already_processed.end())
					continue;
				//Now search of idx_facet in the this list
				for(int p=0;p<n_faces_surrounding;p++)
				{
					int face_id_search = ver.GetIncFacet(p);
					if(face_id_search==idx_facet)
					{
						seed = idx_facet;
						next_seed_found = true;
						break; //p loop breaks
					}
				}
				if(next_seed_found)
					break;//k loop
			}
			// if(!next_seed_found)
			// 	std::cout << "ERRRR";
		}
		// std::cout << "\n";
		(*output)->AddFacet(temp_face_buffer);
	}
}

glm::vec3 catmull_clark_fv(int facet_id,Mesh* mesh)
{
	glm::vec3 x(0.0,0.0,0.0);
	TopoFacet current_face = mesh->GetFacet(facet_id);
	int n_vertices = current_face.GetNumberVertices();
	for(int i=0;i<n_vertices;i++)
	{
		int cur_v_idx = current_face.GetVertexInd(i);
		GeomVert cur_v = mesh->GetGeomVertex(cur_v_idx);
		glm::vec3 p(cur_v.GetCo(0),cur_v.GetCo(1),cur_v.GetCo(2));
		x += p;
	}
	float avg_fac = 1.0f/(float)n_vertices ; 
	x = avg_fac*x;
	return x;
}

glm::vec3 catmull_clark_ev(int edge_id,Mesh* mesh)
{
	TopoEdge edge = mesh->GetEdge(edge_id);
	int v_id = edge.GetVertex(0);
	int w_id = edge.GetVertex(1);
	int f1_id = edge.GetIncFacet(0);
	int f2_id = edge.GetIncFacet(1);
	GeomVert v_g = mesh->GetGeomVertex(v_id);
	GeomVert w_g = mesh->GetGeomVertex(w_id);

	glm::vec3 v(v_g.GetCo(0),v_g.GetCo(1),v_g.GetCo(2));
	glm::vec3 w(w_g.GetCo(0),w_g.GetCo(1),w_g.GetCo(2));
	glm::vec3 vf1 = catmull_clark_fv(f1_id,mesh);
	glm::vec3 vf2 = catmull_clark_fv(f2_id,mesh);

	glm::vec3 ev = v+w+vf1+vf2;
	ev = 0.25f * ev;

	return ev;
}

glm::vec3 mid_point_of_edge(int edge_id,Mesh* mesh)
{
	TopoEdge edge = mesh->GetEdge(edge_id);
	int v_id = edge.GetVertex(0);
	int w_id = edge.GetVertex(1);
	GeomVert v_g = mesh->GetGeomVertex(v_id);
	GeomVert w_g = mesh->GetGeomVertex(w_id);
	glm::vec3 v(v_g.GetCo(0),v_g.GetCo(1),v_g.GetCo(2));
	glm::vec3 w(w_g.GetCo(0),w_g.GetCo(1),w_g.GetCo(2));

	glm::vec3 mp = v+w;
	mp = 0.5f * mp;

	return mp;
}

glm::vec3 catmull_clark_vv(int vertex_id,Mesh* mesh)
{
	TopoVert vert = mesh->GetVertex(vertex_id);
	GeomVert gv = mesh->GetGeomVertex(vertex_id);
	glm::vec3 v(gv.GetCo(0),gv.GetCo(1),gv.GetCo(2)); 

	//Calculate the value of Q
	int n_faces = vert.GetNumberIncFacets();
	glm::vec3 Q(0.0,0.0,0.0);
	for(int i=0;i<n_faces;i++)
	{
		int face_id = vert.GetIncFacet(i);
		glm::vec3 p_q = catmull_clark_fv(face_id,mesh);
		Q += p_q;
	}
	float inv_n_faces = 1.0f/(float)n_faces;
	Q = inv_n_faces*Q;

	//Calculate the value of R
	int n_edges = vert.GetNumberIncEdges();
	glm::vec3 R(0.0,0.0,0.0);
	for(int i=0;i<n_edges;i++)
	{
		int edge_id = vert.GetIncEdge(i);
		glm::vec3 p_r = mid_point_of_edge(edge_id,mesh);
		R += p_r;
	} 
	float inv_n_edges = 1.0f/(float)n_edges;
	R = inv_n_edges*R;

	//Calculate the final values
	float nm3 = (float)(n_edges-3); 
	glm::vec3 res = (inv_n_edges)*Q + (2.0f*inv_n_edges)*R + (nm3*inv_n_edges)*v;

	return res;
}

void catmull_clark(Mesh* input,Mesh** output)
{
	//Algorith 
	//For each vertex
	//Fior each face 
	//Find the two edges it shares 
	//Connect them to make a face 
	int n_vertices = input->GetNumberVertices();
	for(int i=0;i<n_vertices;i++)
	{
		TopoVert cur_vert = input->GetVertex(i);
		glm::vec3 vv = catmull_clark_vv(i,input);
		int n_faces_curvertex = cur_vert.GetNumberIncFacets();
		int n_edges_curvertex = cur_vert.GetNumberIncEdges();
		for(int j=0;j<n_faces_curvertex;j++)
		{
			int fac_idx = cur_vert.GetIncFacet(j);
			glm::vec3 fv = catmull_clark_fv(fac_idx,input);
			//First get the current face being processed
			TopoFacet cur_facet = input->GetFacet(fac_idx);
			int n_edges_in_face = cur_facet.GetNumberEdges();
			//Next put all its edges in a set
			std::set<int> all_edges_in_face;
			for(int k=0;k<n_edges_in_face;k++)
			{
				int edg_idx = cur_facet.GetIncEdge(k);
				all_edges_in_face.insert(edg_idx);
			}
			//Now iterate throught the incident edges of this face and get the two incident ones
			std::set<int>::iterator it;
			int e1_idx,e2_idx;
			bool found_e1=false,found_e2=false;
			for(int k=0;k<n_edges_curvertex;k++)
			{
				int edg_idx = cur_vert.GetIncEdge(k);
				it = all_edges_in_face.find(edg_idx);
				if(it!=all_edges_in_face.end())
				{
					if(!found_e1)
					{
						e1_idx = edg_idx;
						found_e1 = true;
					}
					else if(!found_e2)
					{
						e2_idx = edg_idx;
						found_e2 = true;
					}
				}
				if(found_e2)
					break; //inner k loop for search ends here
			}

			glm::vec3 ev1 = catmull_clark_ev(e1_idx,input);
			glm::vec3 ev2 = catmull_clark_ev(e2_idx,input);

			std::vector<GeomVert> temp_face_buffer;
			temp_face_buffer.push_back(GeomVert(fv.x,fv.y,fv.z));
			temp_face_buffer.push_back(GeomVert(ev1.x,ev1.y,ev1.z));
			temp_face_buffer.push_back(GeomVert(vv.x,vv.y,vv.z));
			temp_face_buffer.push_back(GeomVert(ev2.x,ev2.y,ev2.z));

			(*output)->AddFacet(temp_face_buffer);

		}
	}
}

int third_vertex_detection(int tri_id,int v1,int v2,Mesh* mesh)
{
	TopoFacet tri = mesh->GetFacet(tri_id);
	int v1_tri = tri.GetVertexInd(0);
	int v2_tri = tri.GetVertexInd(1);
	int v3_tri = tri.GetVertexInd(2);
	if(v1_tri==v1 && v2_tri==v2)
		return v3_tri;
	else if(v1_tri==v2 && v2_tri==v1)
		return v3_tri;
	else if(v1_tri==v1 && v3_tri==v2)
		return v2_tri;
	else if(v1_tri==v2 && v3_tri==v1)
		return v2_tri;
	else 
		return v1_tri;
}

glm::vec3 loop_ev(int edge_id,Mesh* mesh)
{
	TopoEdge edge = mesh->GetEdge(edge_id);
	int r_id = edge.GetVertex(0);
	int s_id = edge.GetVertex(1);
	int tri_prs_idx = edge.GetIncFacet(0);
	int tri_qrs_idx = edge.GetIncFacet(1);
	//Detect p_id
	int p_id = third_vertex_detection(tri_prs_idx,r_id,s_id,mesh);
	//Detect q_id
	int q_id = third_vertex_detection(tri_qrs_idx,r_id,s_id,mesh);
	GeomVert p_g = mesh->GetGeomVertex(p_id);
	GeomVert q_g = mesh->GetGeomVertex(q_id);
	GeomVert r_g = mesh->GetGeomVertex(r_id);
	GeomVert s_g = mesh->GetGeomVertex(s_id);

	glm::vec3 p(p_g.GetCo(0),p_g.GetCo(1),p_g.GetCo(2));	
	glm::vec3 q(q_g.GetCo(0),q_g.GetCo(1),q_g.GetCo(2));	
	glm::vec3 r(r_g.GetCo(0),r_g.GetCo(1),r_g.GetCo(2));	
	glm::vec3 s(s_g.GetCo(0),s_g.GetCo(1),s_g.GetCo(2));	

	glm::vec3 res = (0.125f*p) + (0.375f*r) + (0.375f*s) + (0.125f*q);

	return res; 
}

glm::vec3 loop_vv(int vertex_id,Mesh* mesh)
{
	TopoVert ver = mesh->GetVertex(vertex_id);
	int n_adj_edj = ver.GetNumberIncEdges();
	glm::vec3 sum_pi(0.0,0.0,0.0);
	for(int i=0;i<n_adj_edj;i++)
	{
		int edg_idx = ver.GetIncEdge(i);
		TopoEdge edge = mesh->GetEdge(edg_idx);
		int v1_idx = edge.GetVertex(0);
		int v2_idx = edge.GetVertex(1);
		int endpoint_idx;
		if(v1_idx==vertex_id)
			endpoint_idx = v2_idx;
		else
			endpoint_idx = v1_idx;
		GeomVert endpoint_g = mesh->GetGeomVertex(endpoint_idx);
		glm::vec3 pi(endpoint_g.GetCo(0),endpoint_g.GetCo(1),endpoint_g.GetCo(2));
		sum_pi += pi;
	}
	GeomVert v_g = mesh->GetGeomVertex(vertex_id);
	glm::vec3 v(v_g.GetCo(0),v_g.GetCo(1),v_g.GetCo(2));
	float inv_n = 1.0f/(float)n_adj_edj;
	glm::vec3 res = (((1.0f-ALPHA)*inv_n)*sum_pi) + (ALPHA*v);
	return res;
}

void loop(Mesh* input,Mesh** output)
{
	//Loop over all the faces
	//Detecting the stuff correctly
	//Then joining them correctly 
	int n_faces = input->GetNumberFacets();
	for(int i=0;i<n_faces;i++)
	{
		TopoFacet face = input->GetFacet(i);
		//Detect the points and edges correctly
		int p_id = face.GetVertexInd(0);
		int q_id = face.GetVertexInd(1);
		int r_id = face.GetVertexInd(2);
		//Now we must correctly allocate the edges as pq,pr,rp
		int pq_idx=-1,qr_idx=-1,rp_idx=-1;
		int e1_idx = face.GetIncEdge(0);
		TopoEdge e1 = input->GetEdge(e1_idx);
		int e2_idx = face.GetIncEdge(1);
		TopoEdge e2 = input->GetEdge(e2_idx);
		int e3_idx = face.GetIncEdge(2);
		TopoEdge e3 = input->GetEdge(e3_idx);
		//Now we guess whats e1_idx
		int ret1 = third_vertex_detection(i,e1.GetVertex(0),e1.GetVertex(1),input);
		if(ret1==p_id)
			qr_idx = e1_idx;
		else if(ret1==q_id)
			rp_idx = e1_idx;
		else
			pq_idx = e1_idx;
		//Now we guess whats e2_idx
		int ret2 = third_vertex_detection(i,e2.GetVertex(0),e2.GetVertex(1),input);
		if(ret2==p_id)
			qr_idx = e2_idx;
		else if(ret2==q_id)
			rp_idx = e2_idx;
		else
			pq_idx = e2_idx;
		//Now we guess whats e3_idx
		int ret3 = third_vertex_detection(i,e3.GetVertex(0),e3.GetVertex(1),input);
		if(ret3==p_id)
			qr_idx = e3_idx;
		else if(ret3==q_id)
			rp_idx = e3_idx;
		else
			pq_idx = e3_idx;
		if(pq_idx==-1 || rp_idx==-1 || qr_idx==-1)
			std::cout << "[ALERT] Loop Subdivision algorithm is malfunctioning.\n";

		glm::vec3 v_p = loop_vv(p_id,input);
		// GeomVert v_p_g(v_p.x,v_p.y,v_p.z);
		glm::vec3 v_q = loop_vv(q_id,input);
		// GeomVert v_q_g(v_q.x,v_q.y,v_q.z);
		glm::vec3 v_r = loop_vv(r_id,input);
		// GeomVert v_r_g(v_r.x,v_r.y,v_r.z);
		glm::vec3 v_pq = loop_ev(pq_idx,input);
		// GeomVert v_pq_g(v_pq.x,v_pq.y,v_pq.z);
		glm::vec3 v_qr = loop_ev(qr_idx,input);
		// GeomVert v_qr_g(v_qr.x,v_qr.y,v_qr.z);
		glm::vec3 v_rp = loop_ev(rp_idx,input);
		// GeomVert v_rp_g(v_rp.x,v_rp.y,v_rp.z);

		(*output)->AddFacet(v_q.x,v_q.y,v_q.z,v_pq.x,v_pq.y,v_pq.z,v_qr.x,v_qr.y,v_qr.z); //1
		(*output)->AddFacet(v_r.x,v_r.y,v_r.z,v_rp.x,v_rp.y,v_rp.z,v_qr.x,v_qr.y,v_qr.z); //2
		(*output)->AddFacet(v_p.x,v_p.y,v_p.z,v_pq.x,v_pq.y,v_pq.z,v_rp.x,v_rp.y,v_rp.z); //3
		(*output)->AddFacet(v_rp.x,v_rp.y,v_rp.z,v_pq.x,v_pq.y,v_pq.z,v_qr.x,v_qr.y,v_qr.z); //4

	}
}

#endif