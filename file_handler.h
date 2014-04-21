#ifndef FILE_HANDLER_
#define FILE_HANDLER_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include "mesh.h"

#define MAXPOINTS 600
#define SCALING 3

bool process_off_file(char* location,Mesh** ret_obj)
{
	std::ifstream in_file(location);
	float point_buffer_x[MAXPOINTS];
	float point_buffer_y[MAXPOINTS];
	float point_buffer_z[MAXPOINTS];
	if(!in_file)
	{
		std::cout << "[ERROR] Selected OFF File cannot be opened ! \n";
		return false;
	}
	std::string first_iden;
	std::string line;
	int line_counter = 1;

	int n_faces=0,n_vertices=0;
	int v_counter = 0;
	while(getline(in_file,line))
	{
		//getline(in_file,line);
		std::istringstream sst(line);
		if(line_counter==1)
		{
		}
		else if(line_counter==2)
		{
			std::string n_v_str,n_f_str;
			sst >> n_v_str >> n_f_str ;
			n_vertices = atoi(n_v_str.c_str());
			n_faces = atoi(n_f_str.c_str());
		}
		else if(line_counter <= (2+n_vertices))
		{
			std::string x_s,y_s,z_s;
			sst >> x_s >> y_s >> z_s;
			point_buffer_x[v_counter] = atof(x_s.c_str())*SCALING;
			point_buffer_y[v_counter] = atof(y_s.c_str())*SCALING;
			point_buffer_z[v_counter++]= atof(z_s.c_str())*SCALING;
			// std::cout << "Processed vertex number " <<  v_counter << "\n";
		}
		else
		{
			std::vector<GeomVert> temp_facet;
			std::string f_d_s;
			sst >> f_d_s;
			int face_degree = atoi(f_d_s.c_str());
			while(face_degree)
			{
				std::string f_v_s;
				sst >> f_v_s;
				int face_vertex = atoi(f_v_s.c_str());
				// std::cout << face_vertex << " ";
				temp_facet.push_back(GeomVert(point_buffer_x[face_vertex],point_buffer_y[face_vertex],point_buffer_z[face_vertex]));		
				face_degree--;		
			}
			(*ret_obj)->AddFacet(temp_facet);
			// std::cout << "\nAdded a face with size " << temp_facet.size() << "\n";
		}
		line_counter++;
	}
	in_file.close();
	return true;
}


#endif


