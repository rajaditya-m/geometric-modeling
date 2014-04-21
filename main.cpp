#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <deque>
#include <sstream>
#include <fstream>

#ifdef __APPLE__  // include Mac OS X verions of headers
#  include <GL/glew.h>
#  include <GLUT/glut.h>
#  include <OpenGL/gl.h>
#  include <GLUI/glui.h>
#else // non-Mac OS X operating systems
#  include <GL/glew.h>
#  include <GL/glui.h>
#  include <GL/freeglut.h>
#  include <GL/freeglut_ext.h>
#endif  // __APPLE__

#include <glm/gtc/swizzle.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/fast_square_root.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glm/gtx/string_cast.hpp" //Delete later

#include <armadillo>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include "mesh.h"
#include "file_handler.h"
#include "subdiv_algo.h"

#define TOLERANCE 1.0e-5
#define TIMEOUTCNTR 100

int face_render_cntr = 0;

//*******************************************************************
//*This contains all the simplified typedefs we need to call CGAL  *
//*******************************************************************
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_2<Kernel> Point_2;
typedef CGAL::Vector_2<Kernel> Vector_2;
typedef CGAL::Point_3<Kernel> Point_3;
typedef CGAL::Vector_3<Kernel> Vector_3;
typedef CGAL::Segment_2<Kernel> Edge_Segment_2;
typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
typedef Delaunay::Face_iterator face_iter;
typedef Delaunay::Edge_iterator edge_iter;
typedef Delaunay::Vertex_iterator vertex_iter;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;
typedef Traits_2::X_monotone_curve_2 Segment_2;

//*******************************************************************
//*This contains all the information required for the manulpulation *
//*of the screen who knows we maynot even need this one 			*
//*******************************************************************

typedef enum {XFORM_NONE,XFORM_ROTATE,XFORM_SCALE} TransformType;
typedef enum {NONE,CLICK_MODE,SHIFT_MODE} ClickType;
ClickType clicked_type = NONE;
TransformType xform_mode = XFORM_NONE;
GLfloat x_angle=0.0,y_angle=0.0,z_angle=0.0;
int press_x, press_y;
int release_x, release_y;
float scale_size = 1.0;
int clicked_pnt_idx = -1;
int global_mouse_x,global_mouse_y;

//*******************************************************************
//*All the stuff required for GLUI Interface 						*
//*Radio Buttons, Checkbox, Interface 								*
//*******************************************************************
GLUI *glui_window;
int enable_grid_lines = 0;
int enable_control_polygon = 0;
int curve_type = 0;
int twodim_axis = 0;
int local_control_enabled = 0;
int num_subdivisions = 10;
int subdiv_iters = 4;
int nurb_k = 3;
int arbbsp_k = 3;
int junction_point_bc = 0;
int use_x_axis = 0;
int num_slices = 10;
int num_slices_soe = 10;
int num_slices_sweep = 5;
float depth = 6.0;
float x_dir=0.0,y_dir=0.0,z_dir=1.0;
int is_compatible = 0;
int surface_type = 0;
int enable_control_polyhedron = 0;
int subdiv_type = 0;
int enforce_g1_cont = 0;
float arclen_u = 0;
float arclen_w = 0;
float start_u = 0;
float start_w = 0;
float len_u = 0;
float len_w = 0;
int draw_closed_bsp_surface_m = 0;
int draw_closed_bsp_surface_n = 0;
int k_s = 4;
int l_s = 4;
int gen_random_weight_tnurbs = false;
int io_type = 0;
int recons_type = 0;
int always_open_recons = 0;

//*******************************************************************
//*This is the data for control points generated	  				*
//*we stroe them all and then render then							*
//*******************************************************************

std::vector<glm::vec2> control_point_data_xy;
std::vector<glm::vec2> control_point_data_yz;
std::vector<glm::vec2> repos_points_old;

glm::vec2 bezier_curve_data[1000];
glm::vec2 old_bezier_curve_data[1000];
glm::vec2 rational_bezier_curve_data[11];
glm::vec2 cubic_bspline_data[2000];
std::deque<glm::vec2> subdivision_curve_data;
std::vector<glm::vec2> nurbs_data;
std::vector<glm::vec2> arbbsp_data;

//This is a new data strcture made exclusively for the solid generation


int num_points_bezier_curve=0;
int num_points_old_bezier_curve=0;
int num_points_rational_bezier_curve=0;
int num_points_cubicbspline=0;

float weight_rational_bezier[15] ;
int number_of_weights = 0;

float weight_nurbs[15];
int number_of_weights_nurbs=0;

float arbbsp_knot_vector[15];
int number_of_knot_vectors=0;

bool shift_mode_activated=false;
bool add_disp_c1_cont = false;
bool nurbs_ready_for_rendering = false;
bool arbbsp_ready_for_rendering = false;


//We dont need these possibly
int num_control_point = 0; //@TODO Delete

float flat_mat_spline_uniform[16] = {-1.0,3.0,-3.0,1.0,
									3.0,-6.0,0.0,4.0,
									-3.0,3.0,3.0,1.0,
									1.0,0.0,0.0,0.0};
float flat_mat_cubbez[16] = {
	1.0,-3.0,3.0,-1.0,
	0.0,3.0,-6.0,3.0,
	0.0,0.0,3.0,-3.0,
	0.0,0.0,0.0,1.0
};
arma::fmat::fixed<4,4> mat_spline_uniform(flat_mat_spline_uniform);
arma::fmat::fixed<10,4> u_mat_all;
arma::fmat::fixed<4,10> w_mat_all;
arma::fmat::fixed<10,4> cbs_multipler;
arma::fmat::fixed<4,10> cbs_multipler_w;

arma::fmat::fixed<4,4> mat_cubbez(flat_mat_cubbez);
arma::fmat::fixed<9,4> cubbez_multiplier;

bool view_3d = false;

//*******************************************************************
//*These are the data structures for the surface of revolution		*
//*******************************************************************

std::vector<glm::vec3> sor_point_data;
std::vector<int> sor_indices;
int num_quads_sor=0;
bool sor_display  = false;

//*******************************************************************
//*TheQUADSe data structures for the surface of extrusion 		*
//*******************************************************************

std::vector<glm::vec3> soe_point_data;
std::vector<int> soe_indices;
int num_quads_soe=0;
bool soe_display  = false;

std::vector<glm::vec2> history_curve;

//*******************************************************************
//*This is the 3 global mesh daat structures that we shall be using	*
//*******************************************************************
Mesh sor_mesh;
Mesh ext_mesh;

//*******************************************************************
//*These are the data structures for the sweep operations   		*
//*******************************************************************

std::vector<glm::vec3> sweep_point_data;
std::vector<int> sweep_indices;
int num_quads_sweep=0;
bool sweep_display  = false;

//*******************************************************************
//*These are the data structures for the lofting operations   		*
//*******************************************************************

bool lofting_context_mode = false;
std::vector<glm::vec2> saved_curve_samples;
std::vector<glm::vec2> saved_control_points;
int saved_curve_type;
bool loft_display = false;
int num_quads_loft = 0;
std::vector<glm::vec3> loft_point_data;
std::vector<int> loft_indices;

//*******************************************************************
//*These are the data structures for the surface gen operations   	*
//*******************************************************************

int m; //Number of rows
int n; //Number of cols
glm::vec3 control_net_array[100][100];
bool possibly_closed = false;
//*******************************************************************
//*These are the data structures for the beziersurface operations 	*
//*******************************************************************
int num_samples_in_u = 11;
int num_samples_in_v = 11;
bool bez_surface_rendered = false;
glm::vec3 bezier_surface_samples[50][50];

//*******************************************************************
//*These are the data structures for the cubic b-spline operations 	*
//*******************************************************************
int num_samples_in_u_bsp = 0;
int num_samples_in_v_bsp = 0;
bool bsp_surface_rendered = false;
glm::vec3 bspline_surface_samples[500][500];

//*******************************************************************
//*These are the data structures for the NURBS 			operations 	*
//*******************************************************************
int number_of_knot_vectors_i = 0;
int number_of_knot_vectors_j = 0;
float knots_i[30];
float knots_j[30];
float control_weights[100][100];
glm::vec3 tnurbs_surface_samples[500][500];
int tnurbs_samples_in_x=0;
int tnurbs_samples_in_y=0;
bool tnurbs_surface_rendered = false;


//*******************************************************************
//*These are the data structures for the storage of OFF FILE 	*
//*******************************************************************
bool off_exists = false;
Mesh* shape ;
std::vector<glm::vec3> shape_facet_colors;

//*******************************************************************
//*These are the data structures for the subdivions                 *
//*******************************************************************
Mesh* subdiv_res_ds;
Mesh* subdiv_res_cc;
Mesh* subdiv_res_loop;
bool ds_exists = false;
bool cc_exists = false;
bool loop_exists = false;
std::vector<glm::vec3> ds_facet_colors;
std::vector<glm::vec3> cc_facet_colors;
std::vector<glm::vec3> loop_facet_colors;

//*******************************************************************
//*These are the data structures for the recons                 *
//*******************************************************************
std::vector<std::pair <glm::vec2,glm::vec2> > render_recons;

//*******************************************************************
//*This is the data for banner message generated	  				*
//*we stroe them all and then render then							*
//*******************************************************************

std::string banner_array[36] = {
	"Enter more control points to match given number    ", //0
	"Cubic B-Spline Curve is successfully  rendered     ", //1
	"Bezier Curve is successfully rendered              ", //2
	"Mouse Click outside grid - will not be registered  ", //3
	"More Control Points are needed for Cubic B Splines ", //4
	"Excess control point will be ignored               ", //5
	"Bezier curve only allows pseudo local control      ", //6
	"Rational Bezier Curve is successfully rendered     ", //7
	"Weight File successfully read                      ", //8
	"Insufficient weight data to render Rational Bezier ", //9
	"Subdivision curve rendered successfully            ", //10
	"Selected point was successfully duplicated         ", //11
	"No point was selected - duplication not successful ", //12
	"C-1 Continuity Split was successful                ", //13
	"Too few points for C-1 Split (Min : 7)             ", //14
	"NURBS Weight file successfully read                ", //15
	"Insufficient weight data to render NURBS Curve     ", //16
	"NURBS Curve successfully rendered                  ", //17
	"Knot Vector File successfully read                 ", //18
	"Insufficient number of knot vector to render B-Sp. ", //19
	"Arbitrary Degree B-Spline successfully rendered    ", //20
	"2-D Axis changed to Y-Z Axis                       ", //21
	"2-D Axis changed to X-Y Axis                       ", //22
	"Self Intersection detected. No action will happen  ", //23
	"Please draw the next lofting contour               ", //24
	"Bezier surface patch successfully generated        ", //25
	"B-Spline surface successfully generated            ", //26
	"OFF File read successfully to memory.              ", //27
	"Error reading OFF File to memory.                  ", //28
	"No OFF File in memory for selected operation.      ", //29
	"Doo-Sabine Subdivision successful.                 ", //30
	"Catmull-Clark Subdivision successful.              ", //31
	"Loop Subdivision successful.                       ", //32
	"Loop Subdivision requires triangular mesh only.    ", //33
	"Reparameterization successfull.                    ", //34
	"Selected point was successfully deleted.           " //35
};

int banner_warning_type[36] = {
	0,2,2,1,1,2,1,2,2,0,2,2,0,2,0,2,0,2,2,0,2,2,2,0,2,2,2,2,0,0,2,2,2,0,2,2
};

int cur_banner_id = 0;
long long int coefficient_table[15][15];



inline long long int C(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++) {
        ans *= n - r + i;
        ans /= i;
    }
    return ans;
}

float* generate_nonuniform_knot_vector(int k,int n)
{
	//A total of k+n+1 knot vectors will be generated
	float* knot_vector_array = new float[k+n+1];
	for(int i=0;i<k+n+1;i++)
	{
		if(i<k)
			knot_vector_array[i] = 0.0;
		else if(i>n)
			knot_vector_array[i] = (float)(n-k+2);
		else
			knot_vector_array[i] = (float)(i-k+1);
	}

	for(int i=0;i<k+n+1;i++)
	{
		knot_vector_array[i] = knot_vector_array[i]/(float)(n-k+2);
	}

	/*
	for(int i=0;i<k+n+1;i++)
	{
		std::cout << knot_vector_array[i] << " " ;
	}
	std::cout << "\n";
	*/
	return knot_vector_array;
}

float* generate_coefficients_for_u(int k,int n, float* knot_vector,float u)
{
	//First determine the so-called switches for the value of u
	//i.e owhich range it lies in
	float* coeff_vector = new float[n+1];
	float* crawl_space = new float[n+1];
	int nonzero_i=-1;
	for(int i=0;i<k+n;i++)
	{
		if((u>=knot_vector[i]) && (u<knot_vector[i+1]))
		{
			nonzero_i = i;
			break;
		}
	}
	// std::cout << "Non zero detected as :" << nonzero_i << "\n";
	if(nonzero_i==-1)
	{
		//Test if the value of u is 1
		if(fabs(u-1.0)<TOLERANCE)
		{
			// std::cout << "One value detected.\n";
			for(int i=0;i<k+n;i++)
			{
				if(fabs(u-knot_vector[i])<TOLERANCE)
				{
					// std::cout << "Fresh val:" << i << " KV:" <<  knot_vector[i] << "Val:" << abs(u-knot_vector[i]) << "\n";
					nonzero_i = i-1;
					break;
				}
			}
		}
		else
		{
			std::cout << "Something went very very wrong here.";
		}
	}
	// special one for k=0
	for(int i=0;i<n+1;i++)
	{
		if(nonzero_i==i)
			coeff_vector[i]=1.0;
		else
			coeff_vector[i]=0.0;
	}
	// for(int a=0;a<=n;a++)
	// {
	// 	std::cout << coeff_vector[a] << " ";
	// }
	// std::cout << "\n";
	//Now we will do the calculations upto n+1
	//for upto k+1
	for(int j=2;j<=k;j++)
	{
		memcpy(crawl_space,coeff_vector,sizeof(float)*(n+1));
		for(int i=0;i<=n;i++)
		{
			//First check if the denoms are OK
			float denom1 = knot_vector[i+j-1]-knot_vector[i];
			float denom2 = knot_vector[i+j]-knot_vector[i+1];
			float final_val = 0.0;
			if(denom1-0.0>1.0e-5)
			{
				final_val += (((u-knot_vector[i])*crawl_space[i])/(denom1));
			}
			if(denom2-0.0>1.0e-5)
			{
				final_val += (((knot_vector[i+j]-u)*crawl_space[i+1])/(denom2));
			}
			coeff_vector[i] = final_val;
		}
	}
	/*
	std::cout << "NZ:" <<  nonzero_i <<" n:" << n << " k:" << k << " u:" << u << "\n";
	//This should have coeff vector
	for(int a=0;a<=n;a++)
	{
		std::cout << coeff_vector[a] << " ";
	}
	std::cout << "\n";
	*/
	return coeff_vector;
}

bool intersecting_segment_test(glm::vec2 &l1_a,glm::vec2 &l1_b,glm::vec2 &l2_a,glm::vec2 &l2_b)
{
	glm::vec2 p = l1_a;
	glm::vec2 r = l1_b - l1_a;
	glm::vec2 q = l2_a;
	glm::vec2 s = l2_b - l2_a;

	glm::vec3 p3(p,0.0);
	glm::vec3 q3(q,0.0);
	glm::vec3 r3(r,0.0);
	glm::vec3 s3(s,0.0);

	glm::vec3 cp_rs = glm::cross(r3,s3);
	glm::vec3 cp_qmps = glm::cross(q3-p3,s3);

	glm::vec3 cp_pmqr = glm::cross(p3-q3,r3);

	float val_cp_rs = cp_rs.z;
	if(fabs(val_cp_rs-0.0)<TOLERANCE)
		return false;
	else
	{
		float t = cp_qmps.z/cp_rs.z;
		float u = -(cp_pmqr.z/cp_rs.z);

		if(t>0.0 && t<1.0 && u>0.0 && u<1.0)
			return true;
		else
			return false;
	}
}

bool check_for_self_intersections(std::vector<glm::vec2> curve)
{
	int n_samples = curve.size();
	for(int i=0;i<n_samples-1;i++)
	{
		//One segment is i and i+1
		for(int j=i+2;j<n_samples-1;j++)
		{
			if(i==j)
				continue;
			if(intersecting_segment_test(curve[i],curve[i+1],curve[j],curve[j+1]))
				return true;
		}
	}
	return false;
}

void init()
{
	// int* kv = generate_nonuniform_knot_vector(3,5);
	// generate_coefficients_for_u(3,5,kv,1.3);
	//Populate the basis function matrix
	for(int i=0;i<15;++i)
	{
		for(int j=0;j<15;++j)
		{
			coefficient_table[i][j] = C(i,j);
		}
	}
	//Populate the u_mat_all
	u_mat_all.ones();
	w_mat_all.ones();
	float u_val = 0.0;
	float gap = 1.0/9.0;
	for(int i=0;i<10;i++)
	{
		std::cout << "Current u val is " << u_val << "\n";
		u_mat_all(i,0) = u_val * u_val * u_val;
		u_mat_all(i,1) = u_val * u_val;
		u_mat_all(i,2) = u_val;

		w_mat_all(0,i) = u_val * u_val * u_val;
		w_mat_all(1,i) = u_val * u_val;
		w_mat_all(2,i) = u_val;

		u_val += gap;
	}
	//Now generate the mega matrix for super fast generation
	float scalar_multipiler = 1.0/6.0;
	cbs_multipler = u_mat_all * mat_spline_uniform;
	cbs_multipler = scalar_multipiler * cbs_multipler;

	arma::fmat::fixed<4,4> mat_spline_uniform_transpose = mat_spline_uniform.t();


	cbs_multipler_w = mat_spline_uniform_transpose * w_mat_all;
	cbs_multipler_w = scalar_multipiler * cbs_multipler_w;

	arma::fmat::fixed<9,4> t_mat_all;
	t_mat_all.ones();
	u_val = 0.1;
	for(int i=0;i<9;i++)
	{
		t_mat_all(i,3) = u_val * u_val * u_val;
		t_mat_all(i,2) = u_val * u_val;
		t_mat_all(i,1) = u_val;
		u_val += 0.1;
	}
	cubbez_multiplier = t_mat_all * mat_cubbez;
}

glm::vec2 decasteljau_bezier(std::vector<glm::vec2>& control_point,float u)
{
	int n = control_point.size();
	if(n==1)
	{
		return control_point[0];
	}
	else
	{
		std::vector<glm::vec2> new_control_point;
		for(int i=0;i<n-1;i++)
		{
			glm::vec2 temp_point = control_point[i] + (u*(control_point[i+1]-control_point[i]));
			new_control_point.push_back(temp_point);
		}
		return decasteljau_bezier(new_control_point,u);
	}
}

std::deque<glm::vec2> one_subdivide(std::vector<glm::vec2> &point_set,std::deque<glm::vec2> &polygon_1,std::deque<glm::vec2> &polygon_2,float u)
{
	int n = point_set.size();
	if(n==1)
	{
		std::deque<glm::vec2> fin_res(polygon_1);
		fin_res.push_back(point_set[0]);
		fin_res.insert(fin_res.end(),polygon_2.begin(),polygon_2.end());
		return fin_res;
	}
	else
	{
		polygon_1.push_back(point_set[0]);
		polygon_2.push_front(point_set[n-1]);
		std::vector<glm::vec2> new_point_set(n-1);
		for(int i=0;i<=n-2;++i)
		{
			new_point_set[i] = point_set[i] + u*(point_set[i+1]-point_set[i]);
		}
		return one_subdivide(new_point_set,polygon_1,polygon_2,u);
	}
}

std::deque<glm::vec2> subdivide(std::vector<glm::vec2> &point_set,int m,float u)
{
	std::deque<glm::vec2> empty_polygon_1;
	std::deque<glm::vec2> empty_polygon_2;
	if(m==1)
	{
		return one_subdivide(point_set,empty_polygon_1,empty_polygon_2,u);
	}
	else
	{
		std::deque<glm::vec2> ret_val = one_subdivide(point_set,empty_polygon_1,empty_polygon_2,u);
		int ret_val_size = ret_val.size();
		std::deque<glm::vec2>::iterator mid_end = ret_val.begin() + (ret_val_size/2 + 1);
		std::deque<glm::vec2>::iterator mid_start = ret_val.begin() + (ret_val_size/2 );
		std::vector<glm::vec2> point_set_1(ret_val.begin(),mid_end);
		std::vector<glm::vec2> point_set_2(mid_start,ret_val.end());
		std::deque<glm::vec2> ret_val_1 = subdivide(point_set_1,m-1,u);
		std::deque<glm::vec2> ret_val_2 = subdivide(point_set_2,m-1,u);
		for(int i =1;i<ret_val_2.size();i++)
		{
			ret_val_1.push_back(ret_val_2[i]);
		}
		return ret_val_1;
	}
}

std::vector<glm::vec2> get_points_on_cubic_bezier(glm::vec2 &p0,glm::vec2 &p1,glm::vec2 &p2,glm::vec2 &p3)
{
	std::vector<glm::vec2> result(11);
	result[0] = p0;
	arma::fmat::fixed<4,2> P;
	P(0,0) = p0.x; P(0,1) = p0.y;
	P(1,0) = p1.x; P(1,1) = p1.y;
	P(2,0) = p2.x; P(2,1) = p2.y;
	P(3,0) = p3.x; P(3,1) = p3.y;
	arma::fmat::fixed<9,2> R;
	R = cubbez_multiplier * P;
	for(int j=0;j<9;j++)
	{
			glm::vec2 temp_pnt(R(j,0),R(j,1));
			result[j+1] = temp_pnt;
	}
	result[10]=p3;
	return result;
}

void generate_bezier_curve(bool edited=false)
{
	float gap = 1.0/num_subdivisions;
	float u = 0.0;
	int n;
	if(twodim_axis==0)
	 	n = control_point_data_xy.size();
	else
		n = control_point_data_yz.size();

	if(edited && local_control_enabled)
	{
		int adjusted_n = n-1;
		if(clicked_pnt_idx==0)
		{
			if(twodim_axis==0)
				bezier_curve_data[0] = control_point_data_xy[0];
			else
				bezier_curve_data[0] = control_point_data_yz[0];

		}
		else if(clicked_pnt_idx==adjusted_n)
		{
			if(twodim_axis==0)
				bezier_curve_data[num_subdivisions] = control_point_data_xy[adjusted_n];
			else
				bezier_curve_data[num_subdivisions] = control_point_data_yz[adjusted_n];
		}
		else
		{
			u = gap;
  			//Get the value of u which is influenced by this point
			float u_inf = (float)clicked_pnt_idx/(float)adjusted_n;
			//Normalize it within the two u's that will be actually influence
			for(int i=1;i<num_subdivisions;i++)
			{
				float lowlim = u;
				float highlim = u+gap;
				if(u_inf>=lowlim && u_inf<=highlim)
				{
					if(twodim_axis==0)
					{
						bezier_curve_data[i] = decasteljau_bezier(control_point_data_xy,lowlim);
						bezier_curve_data[i+1] = decasteljau_bezier(control_point_data_xy,highlim);
					}
					else
					{
						bezier_curve_data[i] = decasteljau_bezier(control_point_data_yz,lowlim);
						bezier_curve_data[i+1] = decasteljau_bezier(control_point_data_yz,highlim);
					}
					break;
				}
				u += gap;
			}
		}
		cur_banner_id = 6;
		return;
	}
	if(n==0)
		return;
	if(n==1)
	{
		num_points_bezier_curve = 1;
		if(twodim_axis==0)
			bezier_curve_data[0] = control_point_data_xy[0];
		else
			bezier_curve_data[0] = control_point_data_yz[0];
	}
	else if(n==2)
	{
		num_points_bezier_curve = 2;
		if(twodim_axis==0)
		{
			bezier_curve_data[0] = control_point_data_xy[0];
			bezier_curve_data[1] = control_point_data_xy[1];
		}
		else
		{
			bezier_curve_data[0] = control_point_data_yz[0];
			bezier_curve_data[1] = control_point_data_yz[1];
		}
	}
	else
	{
		if(twodim_axis==0)
			bezier_curve_data[0] = control_point_data_xy[0];
		else
			bezier_curve_data[0] = control_point_data_yz[0];
		for(int i=1;i<num_subdivisions;++i)
		{
			u += gap;
			glm::vec2 p_u;
			if(twodim_axis==0)
				p_u = decasteljau_bezier(control_point_data_xy,u);
			else
				p_u = decasteljau_bezier(control_point_data_yz,u);

			bezier_curve_data[i] = p_u;
		}
		if(twodim_axis==0)
		{
			bezier_curve_data[num_subdivisions] = control_point_data_xy[n-1];
		}
		else
		{
			bezier_curve_data[num_subdivisions] = control_point_data_yz[n-1];
		}
		num_points_bezier_curve = num_subdivisions+1;
	}
	cur_banner_id = 2;
}

void generate_rational_bezier_curve()
{
	float gap = 1.0/num_subdivisions;
	float u = 0.0;
	int n ;
	if(twodim_axis==0)
		n = control_point_data_xy.size()-1;
	else
		n = control_point_data_yz.size()-1;
	if(n<1)
		return;
	if(number_of_weights<n+1)
	{
		cur_banner_id = 9;
		num_points_rational_bezier_curve = 0;
		return;
	}
	if(twodim_axis==0)
		rational_bezier_curve_data[0] = control_point_data_xy[0];
	else
		rational_bezier_curve_data[0] = control_point_data_yz[0];

	for(int i=1;i<num_subdivisions;++i)
	{
		u += gap;
		glm::vec2 final_pnt(0.0);
		float den_sum = 0.0;
		for(int j=0;j<=n;++j)
		{
			float exp_val_1 = pow(1.0-u,n-j);
			float exp_val_2 = pow(u,j);
			long long int coeff = coefficient_table[n][j];
			float final_coeff = exp_val_1 * exp_val_2 * (float)coeff * weight_rational_bezier[j];
			glm::vec2 temp;
			if(twodim_axis==0)
				temp = final_coeff * control_point_data_xy[j];
			else
				temp = final_coeff * control_point_data_yz[j];
			final_pnt += temp;
			den_sum += final_coeff;
		}
		float inv_den_sum = 1.0/den_sum;
		rational_bezier_curve_data[i] = inv_den_sum*final_pnt;
	}
	if(twodim_axis==0)
		rational_bezier_curve_data[num_subdivisions] = control_point_data_xy[n];
	else
		rational_bezier_curve_data[num_subdivisions] = control_point_data_yz[n];

	num_points_rational_bezier_curve = num_subdivisions+1;
	cur_banner_id = 7;
}

void generate_cubicbspline_curve(bool edited=false)
{
	int n ;
	if(twodim_axis==0)
		n = control_point_data_xy.size();
	else
		n = control_point_data_yz.size();
	int cbs_counter=0;
	if(n<4)
	{
		cur_banner_id = 4;
		return;
	}
	//Lines of reasoning for the 5 points
	arma::fmat::fixed<4,2> point_mat;
	arma::fmat::fixed<10,2> result_mat;
	if(!edited)
	{
		for(int i=1;i<=n-3;++i)
		{
			if(twodim_axis==0)
			{
				point_mat(0,0) = control_point_data_xy[i-1].x;
				point_mat(0,1) = control_point_data_xy[i-1].y;
				point_mat(1,0) = control_point_data_xy[i].x;
				point_mat(1,1) = control_point_data_xy[i].y;
				point_mat(2,0) = control_point_data_xy[i+1].x;
				point_mat(2,1) = control_point_data_xy[i+1].y;
				point_mat(3,0) = control_point_data_xy[i+2].x;
				point_mat(3,1) = control_point_data_xy[i+2].y;
			}
			else
			{
				point_mat(0,0) = control_point_data_yz[i-1].x;
				point_mat(0,1) = control_point_data_yz[i-1].y;
				point_mat(1,0) = control_point_data_yz[i].x;
				point_mat(1,1) = control_point_data_yz[i].y;
				point_mat(2,0) = control_point_data_yz[i+1].x;
				point_mat(2,1) = control_point_data_yz[i+1].y;
				point_mat(3,0) = control_point_data_yz[i+2].x;
				point_mat(3,1) = control_point_data_yz[i+2].y;
			}
			result_mat = cbs_multipler * point_mat ;
			for(int j=0;j<10;j++)
			{
				glm::vec2 temp_pnt(result_mat(j,0),result_mat(j,1));
				cubic_bspline_data[cbs_counter++] = temp_pnt;
			}
		}
		num_points_cubicbspline = cbs_counter;
	}
	else
	{
		int start_idx;
		// std::cout << "Local control is in effect for point number " << clicked_pnt_idx << "\n";
		//Depending on the clicked point, it will effect 4 segments
		//segment k1 is the point itself
		int k1 = clicked_pnt_idx;
		if(k1>0 && k1 <=n-3)
		{
			if(twodim_axis==0)
			{
				point_mat(0,0) = control_point_data_xy[k1-1].x;
				point_mat(0,1) = control_point_data_xy[k1-1].y;
				point_mat(1,0) = control_point_data_xy[k1].x;
				point_mat(1,1) = control_point_data_xy[k1].y;
				point_mat(2,0) = control_point_data_xy[k1+1].x;
				point_mat(2,1) = control_point_data_xy[k1+1].y;
				point_mat(3,0) = control_point_data_xy[k1+2].x;
				point_mat(3,1) = control_point_data_xy[k1+2].y;
			}
			else
			{
				point_mat(0,0) = control_point_data_yz[k1-1].x;
				point_mat(0,1) = control_point_data_yz[k1-1].y;
				point_mat(1,0) = control_point_data_yz[k1].x;
				point_mat(1,1) = control_point_data_yz[k1].y;
				point_mat(2,0) = control_point_data_yz[k1+1].x;
				point_mat(2,1) = control_point_data_yz[k1+1].y;
				point_mat(3,0) = control_point_data_yz[k1+2].x;
				point_mat(3,1) = control_point_data_yz[k1+2].y;
			}
			result_mat = cbs_multipler * point_mat ;
			start_idx = (k1-1)*10;
			for(int j=0;j<10;j++)
			{
				glm::vec2 temp_pnt(result_mat(j,0),result_mat(j,1));
				cubic_bspline_data[start_idx++] = temp_pnt;
			}
		}
		//segment k2 is one before +1
		int k2 = clicked_pnt_idx+1;
		if(k2>0 && k2 <=n-3)
		{
			if(twodim_axis==0)
			{
				point_mat(0,0) = control_point_data_xy[k2-1].x;
				point_mat(0,1) = control_point_data_xy[k2-1].y;
				point_mat(1,0) = control_point_data_xy[k2].x;
				point_mat(1,1) = control_point_data_xy[k2].y;
				point_mat(2,0) = control_point_data_xy[k2+1].x;
				point_mat(2,1) = control_point_data_xy[k2+1].y;
				point_mat(3,0) = control_point_data_xy[k2+2].x;
				point_mat(3,1) = control_point_data_xy[k2+2].y;
			}
			else
			{
				point_mat(0,0) = control_point_data_yz[k2-1].x;
				point_mat(0,1) = control_point_data_yz[k2-1].y;
				point_mat(1,0) = control_point_data_yz[k2].x;
				point_mat(1,1) = control_point_data_yz[k2].y;
				point_mat(2,0) = control_point_data_yz[k2+1].x;
				point_mat(2,1) = control_point_data_yz[k2+1].y;
				point_mat(3,0) = control_point_data_yz[k2+2].x;
				point_mat(3,1) = control_point_data_yz[k2+2].y;
			}
			result_mat = cbs_multipler * point_mat ;
			start_idx = (k2-1)*10;
			for(int j=0;j<10;j++)
			{
				glm::vec2 temp_pnt(result_mat(j,0),result_mat(j,1));
				cubic_bspline_data[start_idx++] = temp_pnt;
			}
		}
		//segment k3 is one after -1
		int k3 = clicked_pnt_idx-1;
		if(k3>0 && k3 <=n-3)
		{
			if(twodim_axis==0)
			{
				point_mat(0,0) = control_point_data_xy[k3-1].x;
				point_mat(0,1) = control_point_data_xy[k3-1].y;
				point_mat(1,0) = control_point_data_xy[k3].x;
				point_mat(1,1) = control_point_data_xy[k3].y;
				point_mat(2,0) = control_point_data_xy[k3+1].x;
				point_mat(2,1) = control_point_data_xy[k3+1].y;
				point_mat(3,0) = control_point_data_xy[k3+2].x;
				point_mat(3,1) = control_point_data_xy[k3+2].y;
			}
			else
			{
				point_mat(0,0) = control_point_data_yz[k3-1].x;
				point_mat(0,1) = control_point_data_yz[k3-1].y;
				point_mat(1,0) = control_point_data_yz[k3].x;
				point_mat(1,1) = control_point_data_yz[k3].y;
				point_mat(2,0) = control_point_data_yz[k3+1].x;
				point_mat(2,1) = control_point_data_yz[k3+1].y;
				point_mat(3,0) = control_point_data_yz[k3+2].x;
				point_mat(3,1) = control_point_data_yz[k3+2].y;
			}
			result_mat = cbs_multipler * point_mat ;
			start_idx = (k3-1)*10;
			for(int j=0;j<10;j++)
			{
				glm::vec2 temp_pnt(result_mat(j,0),result_mat(j,1));
				cubic_bspline_data[start_idx++] = temp_pnt;
			}
		}
		//segment k4 is two after -2
		int k4 = clicked_pnt_idx-2;
		if(k4>0 && k4 <=n-3)
		{
			if(twodim_axis==0)
			{
				point_mat(0,0) = control_point_data_xy[k4-1].x;
				point_mat(0,1) = control_point_data_xy[k4-1].y;
				point_mat(1,0) = control_point_data_xy[k4].x;
				point_mat(1,1) = control_point_data_xy[k4].y;
				point_mat(2,0) = control_point_data_xy[k4+1].x;
				point_mat(2,1) = control_point_data_xy[k4+1].y;
				point_mat(3,0) = control_point_data_xy[k4+2].x;
				point_mat(3,1) = control_point_data_xy[k4+2].y;
			}
			else
			{
				point_mat(0,0) = control_point_data_yz[k4-1].x;
				point_mat(0,1) = control_point_data_yz[k4-1].y;
				point_mat(1,0) = control_point_data_yz[k4].x;
				point_mat(1,1) = control_point_data_yz[k4].y;
				point_mat(2,0) = control_point_data_yz[k4+1].x;
				point_mat(2,1) = control_point_data_yz[k4+1].y;
				point_mat(3,0) = control_point_data_yz[k4+2].x;
				point_mat(3,1) = control_point_data_yz[k4+2].y;
			}
			result_mat = cbs_multipler * point_mat ;
			start_idx = (k4-1)*10;
			for(int j=0;j<10;j++)
			{
				glm::vec2 temp_pnt(result_mat(j,0),result_mat(j,1));
				cubic_bspline_data[start_idx++] = temp_pnt;
			}
		}
	}
	cur_banner_id = 1;
}

void generate_subdivision_curves()
{
	int n ;
	if(twodim_axis==0)
		n = control_point_data_xy.size();
	else
		n = control_point_data_yz.size();
	if(n>2)
	{
		if(twodim_axis==0)
			subdivision_curve_data = subdivide(control_point_data_xy,subdiv_iters,0.5);
		else
			subdivision_curve_data = subdivide(control_point_data_yz,subdiv_iters,0.5);
		cur_banner_id = 10;
	}
	else if(n==2)
	{
		std::deque<glm::vec2> temp;
		if(twodim_axis==0)
		{
			temp.push_back(control_point_data_xy[0]);
			temp.push_back(control_point_data_xy[1]);
		}
		else
		{
			temp.push_back(control_point_data_yz[0]);
			temp.push_back(control_point_data_yz[1]);
		}
		subdivision_curve_data = temp;
		cur_banner_id = 10;
	}
}

/*Surface Drawing algirithms starts here*/

void generate_soe(int placeholder)
{
	sor_display = false;
	sweep_display = false;
	loft_display = false;
	std::vector<glm::vec2> current_curve_points;
	soe_indices.clear();
	soe_point_data.clear();
	if(curve_type==0)
	{
		for(int i=0;i<num_subdivisions+1;i++)
		{
			current_curve_points.push_back(bezier_curve_data[i]);
		}
	}
	else if(curve_type==1)
	{
		for(int i=0;i<num_points_cubicbspline;i++)
		{
			current_curve_points.push_back(cubic_bspline_data[i]);
		}
	}
	else if(curve_type==2)
	{
		for(int i=0;i<subdivision_curve_data.size();i++)
		{
			current_curve_points.push_back(subdivision_curve_data[i]);
		}
	}
	else if(curve_type==3)
	{
		for(int i=0;i<num_points_rational_bezier_curve;i++)
		{
			current_curve_points.push_back(rational_bezier_curve_data[i]);
		}
	}
	else if(curve_type==4)
	{
		current_curve_points = arbbsp_data;
	}
	else if(curve_type==5)
	{
		current_curve_points = nurbs_data;
	}
	glm::vec3 direction(x_dir,y_dir,z_dir);
	glm::vec3 dir_uvec =  glm::normalize(direction);
	float cur_depth = 0.0;
	num_quads_soe=0;
	int num_points_on_rendered_curve = current_curve_points.size();
	float disp = depth/(float)num_slices_soe;
	for(int i=0;i<num_slices_soe;i++)
	{
		for(int j=0;j<num_points_on_rendered_curve;j++)
		{
			glm::vec3 tp(current_curve_points[j],0.0);
			glm::vec3 np = tp + cur_depth*dir_uvec;
			soe_point_data.push_back(np);
			control_net_array[i][j] = np;
		}
		cur_depth += disp;
	}
	n = num_points_on_rendered_curve;
	m = num_slices_soe;
	std::cout << m << n;
	int counter=0;
	for(int i=0;i<num_slices_soe-1;i++)
	{
		for(int j=0;j<num_points_on_rendered_curve-1;j++)
		{
			int idx1 = i*(num_points_on_rendered_curve)+j;
			int idx2 = (i+1)*(num_points_on_rendered_curve)+j;
			int idx3 = i*(num_points_on_rendered_curve)+(j+1);
			int idx4 = (i+1)*(num_points_on_rendered_curve)+(j+1);
			soe_indices.push_back(idx1);
			soe_indices.push_back(idx3);
			soe_indices.push_back(idx4);
			soe_indices.push_back(idx2);
			num_quads_soe++;
		}
	}
	std::cout << "[INFO] Number of quads in SOE:" << num_quads_soe << "\n";
	std::cout << "[INFO] Number of vertices in SOE:" << soe_point_data.size() << "\n";
	soe_display = true;
	glutPostRedisplay();
}

void generate_sor_solid(int placeholder)
{
	soe_display = false;
	sweep_display = false;
	loft_display = false;
	std::vector<glm::vec2> current_curve_points;
	sor_indices.clear();
	sor_point_data.clear();
	if(curve_type==0)
	{
		for(int i=0;i<num_subdivisions+1;i++)
		{
			current_curve_points.push_back(bezier_curve_data[i]);
		}
	}
	else if(curve_type==1)
	{
		for(int i=0;i<num_points_cubicbspline;i++)
		{
			current_curve_points.push_back(cubic_bspline_data[i]);
		}
	}
	else if(curve_type==2)
	{
		for(int i=0;i<subdivision_curve_data.size();i++)
		{
			current_curve_points.push_back(subdivision_curve_data[i]);
		}
	}
	else if(curve_type==3)
	{
		for(int i=0;i<num_points_rational_bezier_curve;i++)
		{
			current_curve_points.push_back(rational_bezier_curve_data[i]);
		}
	}
	else if(curve_type==4)
	{
		current_curve_points = arbbsp_data;
	}
	else if(curve_type==5)
	{
		current_curve_points = nurbs_data;
	}
	float theta = 360.0/(float)num_slices;
	if(check_for_self_intersections(current_curve_points))
	{
		cur_banner_id = 23;
		glutPostRedisplay();
		sor_display = false;
		// std::cout << "[RERRR] self intersection deteced.!\n";
		return;
	}
	// num_vertices_sor = 0;
	num_quads_sor = 0;
	glm::mat4 rot_mat;
	int num_points_on_rendered_curve = current_curve_points.size();
	std::cout << num_points_on_rendered_curve;
	m = num_points_on_rendered_curve;
	n = num_slices+1;

	for(int i=0;i<num_slices;i++)
	{
		if(use_x_axis)
			rot_mat = glm::rotate(glm::mat4(1.0f),theta*(i+1),glm::vec3(1.0,0.0,0.0));
		else
			rot_mat = glm::rotate(glm::mat4(1.0f),theta*(i+1),glm::vec3(0.0,1.0,0.0));
		for(int j=0;j<num_points_on_rendered_curve;j++)
		{
			glm::vec4 tp(current_curve_points[j],0.0,1.0);
			glm::vec4 rotated_pnt = tp * rot_mat;
			glm::vec3 rp_3 = glm::swizzle<glm::X,glm::Y,glm::Z>(rotated_pnt);
			sor_point_data.push_back(rp_3);
			control_net_array[j][i] = rp_3;
			if(i==0)
			{
				control_net_array[j][n-1] = rp_3;
			}
		}
	}
	possibly_closed = true;
	// std::cout << "[INFO] Number of vertices in SOR : " << num_vertices_sor << "\n";
	int counter=0;
	for(int i=0;i<num_slices;i++)
	{
		for(int j=0;j<num_points_on_rendered_curve-1;j++)
		{
			int idx1 = i*(num_points_on_rendered_curve)+j;
			int idx2 = ((i+1)%num_slices)*(num_points_on_rendered_curve)+j;
			int idx3 = i*(num_points_on_rendered_curve)+(j+1);
			int idx4 = ((i+1)%num_slices)*(num_points_on_rendered_curve)+(j+1);
			sor_indices.push_back(idx1);
			sor_indices.push_back(idx3);
			sor_indices.push_back(idx4);
			sor_indices.push_back(idx2);
			num_quads_sor++;
		}
	}
	std::cout << "[INFO] Number of quads in SOR:" << num_quads_sor << "\n";
	sor_display = true;
	glutPostRedisplay();
}

void generate_loft_surface(int placeholder)
{
	soe_display = false;
	sor_display = false;
	sweep_display = false;

	loft_indices.clear();
	loft_point_data.clear();
	num_quads_loft = 0;

	std::vector<glm::vec2> second_contour_samples;
	if(curve_type==0)
	{
		for(int i=0;i<num_subdivisions+1;i++)
		{
			second_contour_samples.push_back(bezier_curve_data[i]);
		}
	}
	else if(curve_type==1)
	{
		for(int i=0;i<num_points_cubicbspline;i++)
		{
			second_contour_samples.push_back(cubic_bspline_data[i]);
		}
	}
	else if(curve_type==2)
	{
		for(int i=0;i<subdivision_curve_data.size();i++)
		{
			second_contour_samples.push_back(subdivision_curve_data[i]);
		}
	}
	else if(curve_type==3)
	{
		for(int i=0;i<num_points_rational_bezier_curve;i++)
		{
			second_contour_samples.push_back(rational_bezier_curve_data[i]);
		}
	}
	else if(curve_type==4)
	{
		second_contour_samples = arbbsp_data;
	}
	else if(curve_type==5)
	{
		second_contour_samples = nurbs_data;
	}

	//Now we start the correspondence matching
	int num_samples_1 = saved_curve_samples.size();
	int num_samples_2 = second_contour_samples.size();

	std::vector<std::pair<glm::vec3,glm::vec3> > corresponding_values;

	//Assume first curve is at z-axis 0.0
	//Assume that second curve is at z-axis 10.0 (will be changed later) -> this is the one drawn nright now

	//Just do a point to point matching in this case
	if(num_samples_1 == num_samples_2)
	{
		for(int i=0;i<num_samples_1;i++)
		{
			glm::vec3 p1(saved_curve_samples[i],0.0);
			glm::vec3 cand;
			float min_dist_cand = 1000.0;
			for(int j=0;j<num_samples_2;j++)
			{
				glm::vec3 tp2(second_contour_samples[j],10.0);
				float dist = glm::length(tp2-p1);
				if(dist<min_dist_cand)
				{
					min_dist_cand = dist;
					cand = tp2;
				}
			}
			corresponding_values.push_back(std::make_pair(p1,cand));
		}
	}
	else
	{
		//Remember always upsample coz you know thats cool
		if(num_samples_1>num_samples_2)
		{
			for(int i=0;i<num_samples_1;i++)
			{
				glm::vec3 p1(saved_curve_samples[i],0.0);
				glm::vec3 cand;
				float min_dist_cand = 1000.0;
				bool untouched = true;
				for(int j=0;j<num_samples_2-1;j++)
				{
					glm::vec3 tu(second_contour_samples[j],10.0);
					glm::vec3 tv(second_contour_samples[j+1],10.0);
					glm::vec3 pcad;
					glm::vec3 v = tv-tu;
					glm::vec3 w = p1-tu;
					float c1 = glm::dot(w,v);
					float c2 = glm::dot(v,v);
					// if(c1<=0.0)
					// 	pcad = tu;
					// else if(c2 <= c1)
					// 	pcad = tv;
					// else
					// {
						float b = c1/c2;
						pcad = tu + b*v;
					// }
					float dist = glm::length(pcad-p1);
					if(dist<min_dist_cand)
					{
						min_dist_cand = dist;
						cand = pcad;
					}

				}
				corresponding_values.push_back(std::make_pair(p1,cand));
			}
		}
		else
		{
			for(int i=0;i<num_samples_2;i++)
			{
				glm::vec3 p1(second_contour_samples[i],10.0);
				glm::vec3 cand;
				float min_dist_cand = 1000.0;
				bool untouched = true;
				for(int j=0;j<num_samples_1-1;j++)
				{
					glm::vec3 tu(saved_curve_samples[j],0.0);
					glm::vec3 tv(saved_curve_samples[j+1],0.0);
					glm::vec3 pcad;
					glm::vec3 v = tv-tu;
					glm::vec3 w = p1-tu;
					float c1 = glm::dot(w,v);
					float c2 = glm::dot(v,v);
					// if(c1<=0.0)
					// 	pcad = tu;
					// else if(c2 <= c1)
					// 	pcad = tv;
					// else
					// {
						float b = c1/c2;
						pcad = tu + b*v;
					// }
					float dist = glm::length(pcad-p1);
					if(dist<min_dist_cand)
					{
						min_dist_cand = dist;
						cand = pcad;
					}
				}
				corresponding_values.push_back(std::make_pair(cand,p1));
			}
			num_samples_1 = num_samples_2;
		}

	}

	//Now we generate the solid as requested
	float t = 0.0;
	int num_slices_loft = 10;
	float divs = 1.0/(float)num_slices_loft;
	for(int i=0;i<num_slices_loft;i++)
	{
		for(int j=0;j<num_samples_1;j++)
		{
			glm::vec3 origin = corresponding_values[j].first;
			glm::vec3 dir = (corresponding_values[j].second - corresponding_values[j].first);
			glm::vec3 gen_pnt = origin + t*dir;
			loft_point_data.push_back(gen_pnt);
		}
		t += divs;
	}
	//Now we generate the indices vector
	for(int i=0;i<num_slices_loft-1;i++)
	{
		for(int j=0;j<num_samples_1-1;j++)
		{
			int idx1 = i*(num_samples_1)+j;
			int idx2 = (i+1)*(num_samples_1)+j;
			int idx3 = i*(num_samples_1)+(j+1);
			int idx4 = (i+1)*(num_samples_1)+(j+1);
			loft_indices.push_back(idx1);
			loft_indices.push_back(idx3);
			loft_indices.push_back(idx4);
			loft_indices.push_back(idx2);
			num_quads_loft++;
		}
	}
	loft_display = true;
	glutPostRedisplay();
}

void generate_sweep_solid(int placeholder)
{
	soe_display = false;
	sor_display = false;
	loft_display = false;
	//Insert check here to do history curve test
	//And urernt axis check
	sweep_indices.clear();
	sweep_point_data.clear();
	num_quads_sweep=0;
	std::vector<glm::vec2> trajectory_curve ;
	if(curve_type==0)
	{
		for(int i=0;i<num_subdivisions+1;i++)
		{
			trajectory_curve.push_back(bezier_curve_data[i]);
		}
	}
	else if(curve_type==1)
	{
		for(int i=0;i<num_points_cubicbspline;i++)
		{
			trajectory_curve.push_back(cubic_bspline_data[i]);
		}
	}
	else if(curve_type==2)
	{
		for(int i=0;i<subdivision_curve_data.size();i++)
		{
			trajectory_curve.push_back(subdivision_curve_data[i]);
		}
	}
	else if(curve_type==3)
	{
		for(int i=0;i<num_points_rational_bezier_curve;i++)
		{
			trajectory_curve.push_back(rational_bezier_curve_data[i]);
		}
	}
	else if(curve_type==4)
	{
		trajectory_curve = arbbsp_data;
	}
	else if(curve_type==5)
	{
		trajectory_curve = nurbs_data;
	}
	//The generator wire is assumed to lie in the history curve section :)
	int num_points_trajectory = trajectory_curve.size();
	glm::vec3 ip(0.0,0.0,0.0);
	int num_slices_on_gen_curve = 0;
	//Copy the hostory cruve
	std::vector<glm::vec3> generator_curve;
	for(int i=0;i<history_curve.size();i++)
	{
		glm::vec3 tp(history_curve[i],0.0);
		generator_curve.push_back(tp);
	}

	std::vector<glm::vec3> basis_vectors_in_TBN_space;
	glm::vec3 N(1.0,0.0,0.0);
	glm::vec3 tcp1 = glm::vec3(0.0,trajectory_curve[1].y,trajectory_curve[1].x);
	glm::vec3 tcp0 = glm::vec3(0.0,trajectory_curve[0].y,trajectory_curve[0].x);
	glm::vec3 T = glm::normalize(tcp1-tcp0);
	glm::vec3 B = glm::normalize(glm::cross(N,T));

	for(int j=0;j<generator_curve.size();j++)
	{
		glm::vec3 sample_point_3(generator_curve[j].x,generator_curve[j].y,trajectory_curve[0].x);
		glm::vec3 traj_point_3(0.0,trajectory_curve[0].y,trajectory_curve[0].x);

		glm::vec3 sample_in_local = sample_point_3-traj_point_3;

		glm::vec3 tbn_basis;
		tbn_basis.x = glm::dot(sample_in_local,T);
		tbn_basis.y = glm::dot(sample_in_local,B);
		tbn_basis.z = glm::dot(sample_in_local,N);

		basis_vectors_in_TBN_space.push_back(tbn_basis);

		sweep_point_data.push_back(sample_point_3);
	}

	for(int i=1;i<trajectory_curve.size();i++)
	{
		if(i==trajectory_curve.size()-1)
		{
			glm::vec3 tcpn = glm::vec3(0.0,trajectory_curve[i].y,trajectory_curve[i].x);
			glm::vec3 tcpn_1 = glm::vec3(0.0,trajectory_curve[i-1].y,trajectory_curve[i-1].x);
			T = glm::normalize(tcpn-tcpn_1);
		}
		else
		{
			glm::vec3 tcpi = glm::vec3(0.0,trajectory_curve[i].y,trajectory_curve[i].x);
			glm::vec3 tcpim1 = glm::vec3(0.0,trajectory_curve[i-1].y,trajectory_curve[i-1].x);
			glm::vec3 tcpip1 = glm::vec3(0.0,trajectory_curve[i+1].y,trajectory_curve[i+1].x);
			//Weird correction have to implement because of some bozo
			float len_t_left = glm::length(tcpi-tcpim1);
			float len_t_right = glm::length(tcpip1-tcpi);
			if(fabs(len_t_left-0.0)<TOLERANCE)
			{
				T = glm::normalize(tcpip1-tcpi);
			}
			else if(fabs(len_t_right-0.0)<TOLERANCE)
			{
				T = glm::normalize(tcpi-tcpim1);
			}
			else
			{
				glm::vec3 left_t = glm::normalize(tcpi-tcpim1);
				glm::vec3 right_t = glm::normalize(tcpip1-tcpi);
				T = 0.5f*(left_t+right_t);
			}
		}
		B = glm::normalize(glm::cross(N,T));
		for(int j=0;j<generator_curve.size();j++)
		{
			glm::vec3 local_basis_vectors = basis_vectors_in_TBN_space[j];

			glm::vec3 gen_to_trajectory;
			gen_to_trajectory = local_basis_vectors.x*T;
			gen_to_trajectory += local_basis_vectors.y*B;
			gen_to_trajectory += local_basis_vectors.z*N;

			glm::vec3 traj_point_i = glm::vec3(0.0,trajectory_curve[i].y,trajectory_curve[i].x);
			glm::vec3 pnt = traj_point_i + gen_to_trajectory;

			sweep_point_data.push_back(pnt);
		}
		num_slices_on_gen_curve++;
	}
	int num_points_on_rendered_curve = generator_curve.size();
	std::cout << num_slices_on_gen_curve ;
	for(int i=0;i<num_slices_on_gen_curve-1;i++)
	{
		for(int j=0;j<num_points_on_rendered_curve-1;j++)
		{
			int idx1 = i*(num_points_on_rendered_curve)+j;
			int idx2 = (i+1)*(num_points_on_rendered_curve)+j;
			int idx3 = i*(num_points_on_rendered_curve)+(j+1);
			int idx4 = (i+1)*(num_points_on_rendered_curve)+(j+1);
			sweep_indices.push_back(idx1);
			sweep_indices.push_back(idx2);
			sweep_indices.push_back(idx4);
			sweep_indices.push_back(idx3);
			num_quads_sweep++;
		}
	}
	sweep_display = true;
	glutPostRedisplay();
}

void generate_file(int iplaceholder)
{
	std::ofstream file_output("gen.obj");
	if(soe_display)
	{
		if(!is_compatible)
			file_output << soe_point_data.size() << " " << num_quads_soe << "\n";
		for(std::vector<glm::vec3>::iterator it=soe_point_data.begin();it!=soe_point_data.end();it++)
		{
			if(is_compatible)
				file_output << "v " << it->x << " "<< it->y << " " << it->z << "\n";
			else
				file_output << it->x << " "<< it->y << " " << it->z << "\n";
		}
		int every4 = 0;
		for(std::vector<int>::iterator it = soe_indices.begin();it !=soe_indices.end();it++)
		{
			if(every4==0)
			{
				if(is_compatible)
				{
					file_output << "f ";
				}
				else
				{
					file_output << "4 ";
				}
			}

			if(is_compatible)
				file_output << (*it+1) << " ";
			else
				file_output << (*it) << " ";
			every4++;
			if(every4==4)
			{
				file_output << "\n";
				every4 =0;
			}
		}
	}
	else if(sor_display)
	{
		if(!is_compatible)
			file_output << sor_point_data.size() << " " << num_quads_sor << "\n";
		for(std::vector<glm::vec3>::iterator it=sor_point_data.begin();it!=sor_point_data.end();it++)
		{
			if(is_compatible)
				file_output << "v " << it->x << " "<< it->y << " " << it->z << "\n";
			else
				file_output << it->x << " "<< it->y << " " << it->z << "\n";
		}
		int every4 = 0;
		for(std::vector<int>::iterator it = sor_indices.begin();it !=sor_indices.end();it++)
		{
			if(every4==0)
			{
				if(is_compatible)
				{
					file_output << "f ";
				}
				else
				{
					file_output << "4 ";
				}
			}
			if(is_compatible)
				file_output << (*it+1) << " ";
			else
				file_output << (*it) << " ";
			every4++;
			if(every4==4)
			{
				file_output << "\n";
				every4 =0;
			}
		}
	}
	else if(sweep_display)
	{
		if(!is_compatible)
			file_output << sweep_point_data.size() << " " << num_quads_sweep << "\n";
		for(std::vector<glm::vec3>::iterator it=sweep_point_data.begin();it!=sweep_point_data.end();it++)
		{
			if(is_compatible)
				file_output << "v " << it->x << " "<< it->y << " " << it->z << "\n";
			else
				file_output << it->x << " "<< it->y << " " << it->z << "\n";
		}
		int every4 = 0;
		for(std::vector<int>::iterator it = sweep_indices.begin();it !=sweep_indices.end();it++)
		{
			if(every4==0)
			{
				if(is_compatible)
				{
					file_output << "f ";
				}
				else
				{
					file_output << "4 ";
				}
			}
			if(is_compatible)
				file_output << (*it+1) << " ";
			else
				file_output << (*it) << " ";
			every4++;
			if(every4==4)
			{
				file_output << "\n";
				every4 =0;
			}
		}
	}
	file_output.close();
}

glm::vec3 decasteljau_bezier_surface(float u, float w)
{
	//Imagine for now that there is a a global control point array
	//We will fix this restriction later on
	//Also imagine there is a size variable mxn

	std::vector<glm::vec3> p_i_star;
	std::vector<glm::vec3> p_i_star_j_prev;
	std::vector<glm::vec3> p_i_star_j_current;
	int eff_m = m-1;
	int eff_n = n-1;
	for(int i=0;i<=eff_m;i++)
	{
		for(int j=0;j<=eff_n;j++)
		{
			p_i_star_j_prev.push_back(control_net_array[i][j]);
		}
		for(int j=1;j<=eff_n;j++)
		{
			for(int k=0;k<=(eff_n-j);k++)
			{
				glm::vec3 result_point = ((1.0f-w)*p_i_star_j_prev[k]) +
				(w*p_i_star_j_prev[k+1]);
				p_i_star_j_current.push_back(result_point);
			}
			p_i_star_j_prev.clear();
			p_i_star_j_prev = p_i_star_j_current;
			p_i_star_j_current.clear();
		}
		p_i_star.push_back(p_i_star_j_prev[0]);
		p_i_star_j_prev.clear();
	}

	std::vector<glm::vec3> p_i_star_prev;
	std::vector<glm::vec3> p_i_star_current;

	p_i_star_prev = p_i_star;

	for(int j=1;j<=eff_m;j++)
	{
		for(int i=0;i<=(eff_m-j);i++)
		{
			glm::vec3 result_point = ((1.0f-u)*p_i_star_prev[i]) +
			(u*p_i_star_prev[i+1]);
			p_i_star_current.push_back(result_point);
		}
		p_i_star_prev.clear();
		p_i_star_prev = p_i_star_current;
		p_i_star_current.clear();
	}
	return p_i_star_prev[0];
}

void generate_bezier_surface()
{
	if(enforce_g1_cont)
	{

		//Modify three columns
		//0,1,n-2
		for(int i=0;i<m;i++)
		{
			glm::vec3 dir = glm::normalize(control_net_array[i][1]-control_net_array[i][0]);
		  float opp_len = glm::length(control_net_array[i][n-2]-control_net_array[i][0]);
			glm::vec3 new_p = control_net_array[i][0]-(opp_len*dir);
			control_net_array[i][n-2] = new_p;
		}
	}
	float u = 0.0;
	float w = 0.0;
	float u_gap = 1.0f/(float)(num_samples_in_u-1);
	float w_gap = 1.0f/(float)(num_samples_in_v-1);
	for(int i=0;i<num_samples_in_u;i++)
	{
		for(int j=0;j<num_samples_in_v;j++)
		{
			bezier_surface_samples[i][j] = decasteljau_bezier_surface(u,w);
			w += w_gap;
		}
		w = 0.0;
		u += u_gap;
	}
	//Generate the curve length
	float ulen=0.0,vlen=0.0;
	for(int i=0;i<num_samples_in_v-1;i++)
	{
		ulen += glm::length(bezier_surface_samples[0][i]-bezier_surface_samples[0][i+1]);
	}
	arclen_u = ulen;
	for(int i=0;i<num_samples_in_u-1;i++)
	{
		vlen += glm::length(bezier_surface_samples[i][0]-bezier_surface_samples[i+1][0]);
	}
	arclen_w = vlen;
	glui_window->sync_live();
	bez_surface_rendered = true;
}

void generate_bezier_surface(float su,float eu,float sw,float ew)
{
	if(enforce_g1_cont)
	{

		//Modify three columns
		//0,1,n-2
		for(int i=0;i<m;i++)
		{
			glm::vec3 dir = glm::normalize(control_net_array[i][1]-control_net_array[i][0]);
		  float opp_len = glm::length(control_net_array[i][n-2]-control_net_array[i][0]);
			glm::vec3 new_p = control_net_array[i][0]-(opp_len*dir);
			control_net_array[i][n-2] = new_p;
		}
	}
	float u = su;
	float w = sw;
	float u_gap = (eu-su)/(float)(num_samples_in_u-1);
	float w_gap = (ew-sw)/(float)(num_samples_in_v-1);
	for(int i=0;i<num_samples_in_u;i++)
	{
		for(int j=0;j<num_samples_in_v;j++)
		{
			bezier_surface_samples[i][j] = decasteljau_bezier_surface(u,w);
			w += w_gap;
		}
		w = sw;
		u += u_gap;
	}
	//Generate the curve length
	float ulen=0.0,vlen=0.0;
	for(int i=0;i<num_samples_in_v-1;i++)
	{
		ulen += glm::length(bezier_surface_samples[0][i]-bezier_surface_samples[0][i+1]);
	}
	arclen_u = ulen;
	for(int i=0;i<num_samples_in_u-1;i++)
	{
		vlen += glm::length(bezier_surface_samples[i][0]-bezier_surface_samples[i+1][0]);
	}
	arclen_w = vlen;
	glui_window->sync_live();
	bez_surface_rendered = true;
}

void generate_bspline_surface()
{
	int eff_m = m-1;
	int eff_n = n-1;
	int num_iters_to_run_m = eff_m-2;
	if(draw_closed_bsp_surface_m)
		num_iters_to_run_m = eff_m+1;

	int num_iters_to_run_n = eff_n-2;
	if(draw_closed_bsp_surface_n)
		num_iters_to_run_n = eff_n+1;

	for(int s = 1;s<=num_iters_to_run_m;s++)
	{
		for(int t = 1;t<=num_iters_to_run_n;t++)
		{
			arma::fmat::fixed<4,4> point_x;
			arma::fmat::fixed<4,4> point_y;
			arma::fmat::fixed<4,4> point_z;

			int start_idx_x = (s-1)*10;
			int start_idx_y = (t-1)*10;

			int r=0,c=0;

			for(int a=s-1;a<=s+2;a++)
			{
				for(int b=t-1;b<=t+2;b++)
				{
					point_x(r,c) = control_net_array[a%(eff_m+1)][b%(eff_n+1)].x;
					point_y(r,c) = control_net_array[a%(eff_m+1)][b%(eff_n+1)].y;
					point_z(r,c) = control_net_array[a%(eff_m+1)][b%(eff_n+1)].z;
					c++;
				}
				r++;
				c=0;
			}

			arma::fmat::fixed<10,10> result_x;
			arma::fmat::fixed<10,10> result_y;
			arma::fmat::fixed<10,10> result_z;

			result_x = cbs_multipler * point_x * cbs_multipler_w;
			result_y = cbs_multipler * point_y * cbs_multipler_w;
			result_z = cbs_multipler * point_z * cbs_multipler_w;

			//We got a 10x10 +1strip of ar+1rays
			//Now to place them in the correct net location

			r = start_idx_x;
			c = start_idx_y;
			for(int a=0;a<10;a++)
			{
				for(int b=0;b<10;b++)
				{
					glm::vec3 point_temp(result_x(a,b),result_y(a,b),result_z(a,b));
					bspline_surface_samples[r][c] = point_temp;
					c++;
				}
				r++;
				c = start_idx_y;
			}

			num_samples_in_u_bsp = r;
			num_samples_in_v_bsp = c+10;

		}
	}
	//Generate the curve length
	float ulen=0.0,vlen=0.0;
	for(int i=0;i<num_samples_in_v_bsp-1;i++)
	{
		ulen += glm::length(bspline_surface_samples[0][i]-bspline_surface_samples[0][i+1]);
	}
	arclen_u = ulen;
	for(int i=0;i<num_samples_in_u_bsp-1;i++)
	{
		vlen += glm::length(bspline_surface_samples[i][0]-bspline_surface_samples[i+1][0]);
	}
	arclen_w = vlen;
	glui_window->sync_live();
	bsp_surface_rendered = true;
	//std::cout << "Total number of samples generated is " << num_samples_in_u_bsp << " and " << num_samples_in_v_bsp << "\n";
}

void generate_trimmednurbs_surface()
{
	int eff_m = m-1;
	int eff_n = n-1;

	//At this point you will need
	//(m)x(n) weights
	std::ifstream in_file;
	std::string line;
	if(!gen_random_weight_tnurbs)
	{
		in_file.open("tnurbs_weight.dat");
		for(int i=0;i<m;i++)
		{
			getline(in_file,line);
			std::istringstream sst(line);
			for(int j=0;j<n;j++)
			{
				std::string val;
				sst >> val;
				float hval = atof(val.c_str());
				control_weights[i][j] = hval;
			}
			line = "";
		}
		in_file.close();
	}
	else
	{
		for(int i=0;i<m;i++)
		{
			for(int j=0;j<n;j++)
			{
				float hval =  0.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0-0.0)));
				control_weights[i][j] = hval;
			}
		}
	}
	//Now we need some (eff_m+k+1) knot vectors
	int knots_i_req = eff_m + k_s + 1;
	in_file.open("knots_i.dat");
	for(int i=0;i<knots_i_req;i++)
	{
		float kvecval;
		in_file >> kvecval;
		knots_i[i] = kvecval;
	}
	in_file.close();
	//Now we need some (eff_n+l+1) knot vectors
	int knots_j_req = eff_n + l_s + 1;
	in_file.open("knots_j.dat");
	for(int i=0;i<knots_j_req;i++)
	{
		float kvecval;
		in_file >> kvecval;
		knots_j[i] = kvecval;
	}
	in_file.close();

	//Now we will start the main algorithm as needed
	float u = 0.0;
	float w = 0.0;
	float u_gap = 0.05;
	float w_gap = 0.05;
	tnurbs_samples_in_x = tnurbs_samples_in_y = 0;
	while(u<1.05)
	{
		tnurbs_samples_in_y=0;
		float *coeffs_u = generate_coefficients_for_u(k_s,eff_m,knots_i,u);
		while(w<1.05)
		{
			float *coeffs_v = generate_coefficients_for_u(l_s,eff_n,knots_j,w);
			glm::vec3 final_res(0.0,0.0,0.0);
			float denom = 0.0;
			for(int i=0;i<m;i++)
			{
				for(int j=0;j<n;j++)
				{
					final_res += (control_weights[i][j]*coeffs_u[i]*coeffs_v[j])*control_net_array[i][j];
					denom += (control_weights[i][j]*coeffs_u[i]*coeffs_v[j]);
				}
			}
			float denom_inv = 1.0/denom;
			final_res = denom_inv * final_res;
			tnurbs_surface_samples[tnurbs_samples_in_x][tnurbs_samples_in_y] = final_res;
			tnurbs_samples_in_y++;
			w += w_gap;
		}
		u +=u_gap;
		w = 0.0;
		tnurbs_samples_in_x++;
	}

	std::cout  << "NUmber of samples in X " << tnurbs_samples_in_x << "\n";
	std::cout  << "NUmber of samples in Y " << tnurbs_samples_in_y << "\n";
	tnurbs_surface_rendered = true;
	glutPostRedisplay();
}

void generate_face_colors(Mesh* input,std::vector<glm::vec3> &colvec)
{
	colvec.clear();
	int n_facets = input->GetNumberFacets();
	for(int i=0;i<n_facets;i++)
	{
		float r = 0.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0-0.0)));
		float g = 0.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0-0.0)));
		float b = 0.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0-0.0)));
		glm::vec3 col(r,g,b);
		colvec.push_back(col);
	}
}

void clean_circularity(std::vector<Point_2> &ordered_set)
{
	Point_2 first = ordered_set.front();
	Point_2 last = ordered_set.back();
	if(first==last)
	{
		//Delete the last point to clear circularity 
		ordered_set.pop_back();
	}
	/*
	else
	{
		//Find first element 
		std::vector<Point_2>::iterator iter;
		iter = find(ordered_set.begin()+1,ordered_set.end(),first);
		if(iter!=ordered_set.end())
		{

		}
	}
	*/
}

void curve_reconstruction_crust()
{
	//Define the two triangulations
	Delaunay vor_dual;
	Delaunay del_fin;
	//Define a vector of points
	std::vector<Point_2> samples;
	//First construct the delaunay triangulation of the curve
	for(int i=0;i<control_point_data_xy.size();i++)
	{
		glm::vec2 p = control_point_data_xy[i];
		Point_2 cgal_p(p.x,p.y);
		vor_dual.insert(cgal_p);
		del_fin.insert(cgal_p);
		samples.push_back(cgal_p);
	}
	//Now we iterate over all the voronoi vertices
	face_iter fit = vor_dual.faces_begin();
	for(;fit != vor_dual.faces_end(); ++fit)
	{
		Point_2 vor_vertex = vor_dual.dual(fit);
		del_fin.insert(vor_vertex);
	}
	//Now we iterate over all the edges of this delaunay and add those which are in the sample points
  	std::deque<std::pair<Point_2,Point_2> > residues;
  	std::vector<Point_2> ordered_samples;
	edge_iter eit = del_fin.edges_begin();
	int counter = 0;
	bool first_edge = true;
	for(;eit != del_fin.edges_end();++eit)
	{
		Delaunay::Edge e = *eit;
		Delaunay::Vertex_handle v1 = e.first->vertex((e.second + 1) % 3);
		Delaunay::Vertex_handle v2 = e.first->vertex((e.second + 2) % 3);
		// std::cout << "Edge#" << counter << " " << v1->point() << "\n";
		// std::cout << "Edge#" << counter << " " << v2->point() << "\n";
		Point_2 p1 = v1->point();
		Point_2 p2 = v2->point();
		std::vector<Point_2>::iterator it;
		bool match_1_found=false,match_2_found=false;
		for(it=samples.begin();it != samples.end(); ++it)
		{
			Point_2 ref = *it;
			if(ref==p1)
				match_1_found = true;
			else if(ref==p2)
				match_2_found = true;
			if(match_1_found && match_2_found)
				break;
		}
		if(match_1_found && match_2_found)
		{
			if(first_edge)
			{
				ordered_samples.push_back(p1);
				ordered_samples.push_back(p2);
				first_edge = false;
			}
			else
			{
				std::vector<Point_2>::iterator my_iter;
				//Look for point 1
				my_iter = find(ordered_samples.begin(),ordered_samples.end(),p1);
				if(my_iter!=ordered_samples.end())
				{
					//Check if it is the last point
					std::vector<Point_2>::iterator chk_iter = my_iter+1;
					if(chk_iter==ordered_samples.end())
					{
						//It is the last point insert p2
						ordered_samples.push_back(p2);
					}
					else if(my_iter==ordered_samples.begin())
					{
						ordered_samples.insert(ordered_samples.begin(),p2);
					}
				}
				else
				{
					my_iter = find(ordered_samples.begin(),ordered_samples.end(),p2);
					if(my_iter!=ordered_samples.end())
					{
						std::vector<Point_2>::iterator chk_iter = my_iter+1;
						if(chk_iter==ordered_samples.end())
						{
							//It is the last point insert p1
							ordered_samples.push_back(p1);
						}
						else if(my_iter==ordered_samples.begin())
						{
							ordered_samples.insert(ordered_samples.begin(),p1);
						}
					}
					else
					{
						residues.push_back(std::make_pair(p1,p2));
					}
				}
			}
		}
		counter++;
	}
	//Now we handle the residues correctly
	int time_out = 0;
	while(!residues.empty())
	{
		std::pair<Point_2,Point_2> cur_res = residues.front();
		residues.pop_front();
		Point_2 p1 = cur_res.first;
		Point_2 p2 = cur_res.second;
		Point_2 head = ordered_samples.front();
		Point_2 tail = ordered_samples.back();
		if(head==p1)
		{
			ordered_samples.insert(ordered_samples.begin(),p2);
		}
		else if(tail==p1)
		{
			ordered_samples.push_back(p2);
		}
		else if(head==p2)
		{
			ordered_samples.insert(ordered_samples.begin(),p1);
		}
		else if(tail==p2)
		{
			ordered_samples.push_back(p1);
		}
		else
		{
			residues.push_back(cur_res);
		}
		time_out++;
		if(time_out==TIMEOUTCNTR)
			break;
	}
	render_recons.clear();
	for(int j=0;j<residues.size();j++)
	{
		std::pair <Point_2,Point_2> pr = residues[j];
		Point_2 p1 = pr.first;
		Point_2 p2 = pr.second;
		glm::vec2 tmp1(p1.x(),p1.y());
		glm::vec2 tmp2(p2.x(),p2.y());
		render_recons.push_back(std::make_pair(tmp1,tmp2));
	}
	if(always_open_recons)
	{
		clean_circularity(ordered_samples);
	}
	control_point_data_xy.clear();
	std::vector<Point_2>::iterator it;
	for(it=ordered_samples.begin();it!=ordered_samples.end();++it)
	{
		std::cout << "[" << *it << "]" << "\n";
		glm::vec2 tmp(it->x(),it->y());
		control_point_data_xy.push_back(tmp);
	}
}


void curve_reconstruction_crust_v2()
{
	//Define the two triangulations
	Delaunay vor_dual;
	Delaunay del_fin;
	//Define a vector of points
	std::vector<Point_2> samples;
	//Arragements vector
	Arrangement_2 arr;
	Walk_pl pl(arr);
	//First construct the delaunay triangulation of the curve
	for(int i=0;i<control_point_data_xy.size();i++)
	{
		glm::vec2 p = control_point_data_xy[i];
		Point_2 cgal_p(p.x,p.y);
		vor_dual.insert(cgal_p);
		del_fin.insert(cgal_p);
		samples.push_back(cgal_p);
	}
	//Now we iterate over all the voronoi vertices
	face_iter fit = vor_dual.faces_begin();
	for(;fit != vor_dual.faces_end(); ++fit)
	{
		Point_2 vor_vertex = vor_dual.dual(fit);
		del_fin.insert(vor_vertex);
	}
	//Now we iterate over all the edges of this delaunay and add those which are in the sample points
  	std::deque<std::pair<Point_2,Point_2> > residues;
  	std::vector<Point_2> ordered_samples;
	edge_iter eit = del_fin.edges_begin();
	int counter = 0;
	bool first_edge = true;
	for(;eit != del_fin.edges_end();++eit)
	{
		Delaunay::Edge e = *eit;
		Delaunay::Vertex_handle v1 = e.first->vertex((e.second + 1) % 3);
		Delaunay::Vertex_handle v2 = e.first->vertex((e.second + 2) % 3);
		// std::cout << "Edge#" << counter << " " << v1->point() << "\n";
		// std::cout << "Edge#" << counter << " " << v2->point() << "\n";
		Point_2 p1 = v1->point();
		Point_2 p2 = v2->point();
		std::vector<Point_2>::iterator it;
		bool match_1_found=false,match_2_found=false;
		for(it=samples.begin();it != samples.end(); ++it)
		{
			Point_2 ref = *it;
			if(ref==p1)
				match_1_found = true;
			else if(ref==p2)
				match_2_found = true;
			if(match_1_found && match_2_found)
				break;
		}
		if(match_1_found && match_2_found)
		{
			Segment_2 temp_seg(p1,p2);
  		insert_non_intersecting_curve(arr,temp_seg,pl);
		}
	}
	Arrangement_2::Face_const_iterator fiter;
	for (fiter = arr.faces_begin(); fiter != arr.faces_end(); ++fiter)
	{
		Arrangement_2::Face_const_handle f = fiter;
		if(f->is_unbounded())
		{
			std::cout << "Unbounded face determination.\n";
			Arrangement_2::Hole_const_iterator hi;
  		hi = f->holes_begin();
			Arrangement_2::Ccb_halfedge_const_circulator circ = *hi;
			Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
  		// std::cout << "(" << curr->source()->point() << ")";
  		do
  		{
    		std::cout << "   [" << curr->curve() << "] ";//  "      << "(" << curr->target()->point() << ")";
  		} while (++curr != circ);
  		std::cout << std::endl;
		}
		else
		{
			Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
			Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
  		std::cout << "(" << curr->source()->point() << ")";
  		do
  		{
    		std::cout << "   [" << curr->curve() << "]   "      << "(" << curr->target()->point() << ")";
  		} while (++curr != circ);
  		std::cout << std::endl;
		}
	}
 //    // print_face (fit);
	// // fit = arr.faces_begin();
	// control_point_data_xy.clear();
	// std::vector<Point_2>::iterator it;
	// for(it=ordered_samples.begin();it!=ordered_samples.end();++it)
	// {
	// 	std::cout << "[" << *it << "]" << "\n";
	// 	glm::vec2 tmp(it->x(),it->y());
	// 	control_point_data_xy.push_back(tmp);
	// }
}


/*
void curve_reconstruction_NN_crust()
{
	//Define the solution vectors
	std::deque<std::pair<Point_2,Point_2> > residues;
  std::vector<Point_2> ordered_samples;
  bool first_edge = true;
	//Define the delaunay triangulations
	Delaunay del_fin;
	//Define a vector of points
	std::vector<Point_2> samples;
	//First construct the delaunay triangulation of the curve
	for(int i=0;i<control_point_data_xy.size();i++)
	{
		glm::vec2 p = control_point_data_xy[i];
		Point_2 cgal_p(p.x,p.y);
		del_fin.insert(cgal_p);
		samples.push_back(cgal_p);
	}

	//Now we will traverse every vertex of the triangulation
	vertex_iter vit = del_fin.vertices_begin();
	int pc=0;
	for(;vit != del_fin.vertices_end(); ++vit)
	{
		Point_2 p = vit->point();
		//First we determine the shortest edge and call it pq
		Delaunay::Segment shortest_seg;
		double shortest_len = 1000000.0f;
		int ec =0;
		Delaunay::Edge_circulator start = del_fin.incident_edges(vit);
		Delaunay::Edge_circulator cur = start;
		do
		{
			if(del_fin.is_infinite(cur))
			{
				cur++;
				continue;
			}
			Delaunay::Segment seg = del_fin.segment(cur);
			double len_sq = seg.squared_length();
			double len = sqrt(len_sq);
			if(len<shortest_len)
			{
				shortest_seg = seg;
				shortest_len =len;
			}
			++cur;
			ec++;
		}while(cur!=start);
		Point_2 p1 = shortest_seg.source();
		Point_2 p2 = shortest_seg.target();
		Point_2 q;
		if(p1==p)
			q = p2;
		else
			q = p1;
		//Next we determine the shortest edge ps such that angle pqs is more than 90
		shortest_len = 1000000.0f;
		start = del_fin.incident_edges(vit);
		cur = start;
		bool found_candidate = false;
		do
		{
			if(del_fin.is_infinite(cur))
			{
				cur++;
				continue;
			}
			Delaunay::Segment seg = del_fin.segment(cur);
			Point_2 p1_t = seg.source();
			Point_2 p2_t = seg.target();
			Point_2 cand_s;
			if(p1_t==p)
				cand_s = p2_t;
			else
				cand_s = p1_t;
			CGAL::Angle angle_estimate = CGAL::angle(q,p,cand_s); 
			if(angle_estimate!=CGAL::ACUTE)
			{
				found_candidate = true;
				double len_sq = seg.squared_length();
				double len = sqrt(len_sq);
				if(len<shortest_len)
				{
					shortest_seg = seg;
					shortest_len = len;
				}
			}
			++cur;
		}while(cur!=start);
		if(!found_candidate)
			continue;
		p1 = shortest_seg.source();
		p2 = shortest_seg.target();
		Point_2 s;
		if(p1==p)
			s = p2;
		else
			s = p1;

		//Now we will put edges pq,ps in our map as before 
		std::cout << "Reconst edge [" << q << "] [" << p << "] [" << s << "]\n"; 
		if(first_edge)
		{
			ordered_samples.push_back(q);
			ordered_samples.push_back(p);
			ordered_samples.push_back(s);
			first_edge = false;
			std::cout << "Added edges Q[" << q << "] P[" << p << "] S[" << s << "]\n";
		}
		else
		{
			std::vector<Point_2>::iterator my_iter;
			//Insert edge ps
			//Look for point s 
			my_iter = find(ordered_samples.begin(),ordered_samples.end(),s);
			if(my_iter!=ordered_samples.end())
			{
				//Check if it is the last point
				std::vector<Point_2>::iterator next_elem = my_iter+1;
				std::vector<Point_2>::iterator prev_elem = my_iter-1;
				if(next_elem==ordered_samples.end())
				{
					if(*prev_elem!=p)
					{
						ordered_samples.push_back(p);
						std::cout << "Added edges P[" << p << "] S[" << s << "]\n";
					}
				}
				else if(my_iter==ordered_samples.begin())
				{
					if(*next_elem!=p)
					{
						ordered_samples.insert(ordered_samples.begin(),p);
						std::cout << "Added edges P[" << p << "] S[" << s << "]\n";
					}
				}
			}
			else
			{
				my_iter = find(ordered_samples.begin(),ordered_samples.end(),p);
				if(my_iter!=ordered_samples.end())
				{
					std::vector<Point_2>::iterator next_elem = my_iter+1;
					std::vector<Point_2>::iterator prev_elem = my_iter-1;
					if(next_elem==ordered_samples.end())
					{
						if(*prev_elem!=s)
						{
							ordered_samples.push_back(s);
							std::cout << "Added edges P[" << p << "] S[" << s << "]\n";
						}

					}
					else if(my_iter==ordered_samples.begin())
					{
						if(*next_elem!=s)
						{
							ordered_samples.insert(ordered_samples.begin(),s);
							std::cout << "Added edges P[" << p << "] S[" << s << "]\n";
						}
					}
				}
				else
				{
					residues.push_back(std::make_pair(p,s));
					std::cout << "Residuals P[" << p << "] S[" << s << "]\n";
				}
			}
			//Insert edge qp
			//Look for point q 
			my_iter = find(ordered_samples.begin(),ordered_samples.end(),q);
			if(my_iter!=ordered_samples.end())
			{
				//Check if it is the last point
				std::vector<Point_2>::iterator next_elem = my_iter+1;
				std::vector<Point_2>::iterator prev_elem = my_iter-1;
				if(next_elem==ordered_samples.end())
				{
					if(*prev_elem!=p)
					{
						ordered_samples.push_back(p);
						std::cout << "Added edges P[" << p << "] Q[" << q << "]\n";
					}
				}
				else if(my_iter==ordered_samples.begin())
				{
					if(*next_elem!=p)
					{
						ordered_samples.insert(ordered_samples.begin(),p);
						std::cout << "Added edges P[" << p << "] Q[" << q << "]\n";
					}
				}
			}
			else
			{
				my_iter = find(ordered_samples.begin(),ordered_samples.end(),p);
				if(my_iter!=ordered_samples.end())
				{
					std::vector<Point_2>::iterator next_elem = my_iter+1;
					std::vector<Point_2>::iterator prev_elem = my_iter-1;
					if(next_elem==ordered_samples.end())
					{
						if(*prev_elem!=q)
						{
							ordered_samples.push_back(q);
							std::cout << "Added edges P[" << p << "] Q[" << q << "]\n";
						}
					}
					else if(my_iter==ordered_samples.begin())
					{
						if(*next_elem!=q)
						{
							ordered_samples.insert(ordered_samples.begin(),q);
							std::cout << "Added edges P[" << p << "] Q[" << q << "]\n";
						}
					}
				}
				else
				{
					residues.push_back(std::make_pair(p,q));
					std::cout << "Residuals P[" << p << "] Q[" << q << "]\n";
				}
			}
		}


		pc++;
	}
	std::cout <<  "Finished processing all points\n";
	std::cout << "Size of residue " << residues.size() << "\n";
	std::cout << "Size of reconstructed point samples so far :" << ordered_samples.size() << "\n";
	std::cout << "Original size of samples : " << control_point_data_xy.size() << "\n";
	//Now time to process the residues
	int time_out = 0;
	while(!residues.empty())
	{
		std::pair<Point_2,Point_2> cur_res = residues.front();
		residues.pop_front();
		Point_2 p1 = cur_res.first;
		Point_2 p2 = cur_res.second;
		Point_2 head = ordered_samples.front();
		Point_2 tail = ordered_samples.back();
		if(head==p1)
		{
			ordered_samples.insert(ordered_samples.begin(),p2);
		}
		else if(tail==p1)
		{
			ordered_samples.push_back(p2);
		}
		else if(head==p2)
		{
			ordered_samples.insert(ordered_samples.begin(),p1);
		}
		else if(tail==p2)
		{
			ordered_samples.push_back(p1);
		}
		else
		{
			residues.push_back(cur_res);
		}
		time_out++;
		if(time_out==50)
			break;
	}
	
	// std::cout << "Point set generated \n";
	control_point_data_xy.clear();
	std::vector<Point_2>::iterator it;
	for(it=ordered_samples.begin();it!=ordered_samples.end();++it)
	{
		// std::cout << "[" << *it << "]" << "\n";
		glm::vec2 tmp(it->x(),it->y());
		control_point_data_xy.push_back(tmp);
	}	
	// std::cout << "-----------------------\n";
}
*/

bool add_to_edge_soup(std::vector<std::pair<Point_2,Point_2> > &soup,Point_2 p, Point_2 q)
{
	int num_edges = soup.size();
	bool found = false;
	for(int i = 0;i<num_edges;i++)
	{
		std::pair<Point_2,Point_2> cur_pair = soup[i];
		Point_2 p1 = cur_pair.first;
		Point_2 p2 = cur_pair.second;
		if((p1==p && p2==q) || (p1==q && p2==p))
			found = true;
		if(found)
			break;
	}
	if(!found)
		soup.push_back(std::make_pair(p,q));
	return !found;
}

void curve_reconstruction_NN_crust()
{
	//Define the solution vectors
	std::deque<std::pair<Point_2,Point_2> > residues;
  std::vector<Point_2> ordered_samples;
  std::vector<std::pair<Point_2,Point_2> > edge_soup;
  bool first_edge = true;
	//Define the delaunay triangulations
	Delaunay del_fin;
	//Define a vector of points
	std::vector<Point_2> samples;
	//First construct the delaunay triangulation of the curve
	for(int i=0;i<control_point_data_xy.size();i++)
	{
		glm::vec2 p = control_point_data_xy[i];
		Point_2 cgal_p(p.x,p.y);
		del_fin.insert(cgal_p);
		samples.push_back(cgal_p);
	}

	//Now we will traverse every vertex of the triangulation
	vertex_iter vit = del_fin.vertices_begin();
	int pc=0;
	for(;vit != del_fin.vertices_end(); ++vit)
	{
		Point_2 p = vit->point();
		//First we determine the shortest edge and call it pq
		Delaunay::Segment shortest_seg;
		double shortest_len = 1000000.0f;
		int ec =0;
		Delaunay::Edge_circulator start = del_fin.incident_edges(vit);
		Delaunay::Edge_circulator cur = start;
		do
		{
			if(del_fin.is_infinite(cur))
			{
				cur++;
				continue;
			}
			Delaunay::Segment seg = del_fin.segment(cur);
			double len_sq = seg.squared_length();
			double len = sqrt(len_sq);
			if(len<shortest_len)
			{
				shortest_seg = seg;
				shortest_len =len;
			}
			++cur;
			ec++;
		}while(cur!=start);
		Point_2 p1 = shortest_seg.source();
		Point_2 p2 = shortest_seg.target();
		Point_2 q;
		if(p1==p)
			q = p2;
		else
			q = p1;
		//Next we determine the shortest edge ps such that angle pqs is more than 90
		shortest_len = 1000000.0f;
		start = del_fin.incident_edges(vit);
		cur = start;
		bool found_candidate = false;
		do
		{
			if(del_fin.is_infinite(cur))
			{
				cur++;
				continue;
			}
			Delaunay::Segment seg = del_fin.segment(cur);
			Point_2 p1_t = seg.source();
			Point_2 p2_t = seg.target();
			Point_2 cand_s;
			if(p1_t==p)
				cand_s = p2_t;
			else
				cand_s = p1_t;
			CGAL::Angle angle_estimate = CGAL::angle(q,p,cand_s); 
			if(angle_estimate!=CGAL::ACUTE)
			{
				found_candidate = true;
				double len_sq = seg.squared_length();
				double len = sqrt(len_sq);
				if(len<shortest_len)
				{
					shortest_seg = seg;
					shortest_len = len;
				}
			}
			++cur;
		}while(cur!=start);
		if(!found_candidate)
			continue;
		p1 = shortest_seg.source();
		p2 = shortest_seg.target();
		Point_2 s;
		if(p1==p)
			s = p2;
		else
			s = p1;
		bool res1 = add_to_edge_soup(edge_soup,p,q);
		bool res2 = add_to_edge_soup(edge_soup,p,s);
		pc++;
	}
	std::cout <<  "Finished processing all points\n";
	std::cout << "Trying to reconstruct the sample set.\n";
	std::cout << "Edge soup has exactly " <<  edge_soup.size() << " edges\n";
	//Now we process the edge soup and the residues
	for(int i=0;i<edge_soup.size();i++)
	{
		std::pair<Point_2,Point_2> cur_edge = edge_soup[i];
		Point_2 p1 = cur_edge.first;
		Point_2 p2 = cur_edge.second;
		if(first_edge)
		{
			ordered_samples.push_back(p1);
			ordered_samples.push_back(p2);
			first_edge = false;
		}
		else
		{
			std::vector<Point_2>::iterator my_iter;
			//Look for point 1
			my_iter = find(ordered_samples.begin(),ordered_samples.end(),p1);
			if(my_iter!=ordered_samples.end())
			{
				//Check if it is the last point
				std::vector<Point_2>::iterator chk_iter = my_iter+1;
				if(chk_iter==ordered_samples.end())
				{
					//It is the last point insert p2
					ordered_samples.push_back(p2);
				}
				else if(my_iter==ordered_samples.begin())
				{
					ordered_samples.insert(ordered_samples.begin(),p2);
				}
			}
			else
			{
				my_iter = find(ordered_samples.begin(),ordered_samples.end(),p2);
				if(my_iter!=ordered_samples.end())
				{
					std::vector<Point_2>::iterator chk_iter = my_iter+1;
					if(chk_iter==ordered_samples.end())
					{
						//It is the last point insert p1
						ordered_samples.push_back(p1);
					}
					else if(my_iter==ordered_samples.begin())
					{
						ordered_samples.insert(ordered_samples.begin(),p1);
					}
				}
				else
				{
					residues.push_back(std::make_pair(p1,p2));
				}
			}
		}
	}
	std::cout << "Trying to process the " << residues.size() << " residues\n";
	//Now time to process the residues
	int time_out = 0;
	while(!residues.empty())
	{
		std::pair<Point_2,Point_2> cur_res = residues.front();
		residues.pop_front();
		Point_2 p1 = cur_res.first;
		Point_2 p2 = cur_res.second;
		Point_2 head = ordered_samples.front();
		Point_2 tail = ordered_samples.back();
		if(head==p1)
		{
			ordered_samples.insert(ordered_samples.begin(),p2);
		}
		else if(tail==p1)
		{
			ordered_samples.push_back(p2);
		}
		else if(head==p2)
		{
			ordered_samples.insert(ordered_samples.begin(),p1);
		}
		else if(tail==p2)
		{
			ordered_samples.push_back(p1);
		}
		else
		{
			residues.push_back(cur_res);
		}
		time_out++;
		if(time_out==TIMEOUTCNTR)
			break;
	}
	if(always_open_recons)
	{
		clean_circularity(ordered_samples);
	}
	render_recons.clear();
	for(int j=0;j<residues.size();j++)
	{
		std::pair <Point_2,Point_2> pr = residues[j];
		Point_2 p1 = pr.first;
		Point_2 p2 = pr.second;
		glm::vec2 tmp1(p1.x(),p1.y());
		glm::vec2 tmp2(p2.x(),p2.y());
		render_recons.push_back(std::make_pair(tmp1,tmp2));
	}
	// std::cout << "Point set generated \n";
	control_point_data_xy.clear();
	std::vector<Point_2>::iterator it;
	for(it=ordered_samples.begin();it!=ordered_samples.end();++it)
	{
		// std::cout << "[" << *it << "]" << "\n";
		glm::vec2 tmp(it->x(),it->y());
		control_point_data_xy.push_back(tmp);
	}	
	// std::cout << "-----------------------\n";
}


//*******************************************************************
//*This is the routines for all the drawing     	  				*
//*functions that my code uses          							*
//*******************************************************************

void draw_axis_marks(short int axis,glm::vec3 &col,float center_mark,bool enable_rot = false)
{

	float width = 0.3;

	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(4.5);

	glBegin(GL_LINES);
	glColor3f(col.x,col.y,col.z);
	if(axis==1)
	{
		glVertex3f(0.0,center_mark,0.0);
		glVertex3f(width,center_mark,0.0);
	}
	else if(axis==2)
	{
		glVertex3f(center_mark,0.0,0.0);
		glVertex3f(center_mark,width,0.0);
	}
	else if(axis==3)
	{
		glVertex3f(0.0,0.0,center_mark);
		glVertex3f(0.0,width,center_mark);
	}
	glEnd();

	glPopMatrix();
}

void draw_axes(bool enable_rot = false,int axis_id=0)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(4.5);

	glBegin(GL_LINES);

	if(axis_id==0)
		glColor3f(1.0,0.0,0.0); //red
	else if(axis_id==1)
		glColor3f(0.0,1.0,0.0); //green
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(20.0,0,0);

	glColor3f(0.0,0.0,1.0); //blue
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(0.0,20.0,0);

	if(axis_id==0)
		glColor3f(0.0,1.0,0.0);
	else
		glColor3f(1.0,0.0,0.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(0.0,0,2.0);

	glEnd();

	glPopMatrix();

	glm::vec3 red_color(1.0,0.0,0.0);
	glm::vec3 blue_color(0.0,0.0,1.0);
	glm::vec3 green_color(0.0,1.0,0.0);
	for(int i=1;i<21;i++)
	{

		if(axis_id==0)
			draw_axis_marks(2,red_color,(float)i,enable_rot);
		else
			draw_axis_marks(2,green_color,(float)i,enable_rot);

		draw_axis_marks(1,blue_color,(float)i,enable_rot);
	}
}

void draw_axes_3d(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}
	glLineWidth(4.5);

	glBegin(GL_LINES);

	glColor3f(1.0,0.0,0.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(20.0,0,0);

	glColor3f(0.0,0.0,1.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(0.0,20.0,0);

	glColor3f(0.0,1.0,0.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(0.0,0,20.0);

	glEnd();

	glPopMatrix();

	glm::vec3 red_color(1.0,0.0,0.0);
	glm::vec3 blue_color(0.0,0.0,1.0);
	glm::vec3 green_color(0.0,1.0,0.0);
	for(int i=1;i<21;i++)
	{
		draw_axis_marks(1,blue_color,(float)i,enable_rot);
		draw_axis_marks(2,red_color,(float)i,enable_rot);
		draw_axis_marks(3,green_color,(float)i,enable_rot);
	}
}

void draw_control_polygon(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glLineStipple(1,0x3F07);
	glEnable(GL_LINE_STIPPLE);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,1.0);

	for(int p=0;p<control_point_data_xy.size();++p)
	{
		glVertex3f(control_point_data_xy[p].x,control_point_data_xy[p].y,0.0);
	}

	glEnd();

	glDisable(GL_LINE_STIPPLE);

	glPopMatrix();
}

void draw_residues(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glLineStipple(1,0x3F07);
	glEnable(GL_LINE_STIPPLE);

	glColor3f(0.0,1.0,1.0);

	for(int i=0;i<render_recons.size();i++)
	{
		glm::vec2 p1 = render_recons[i].first;
		glm::vec2 p2 = render_recons[i].second;
		glBegin(GL_LINES);
		glVertex3f(p1.x,p1.y,0.0);
		glVertex3f(p2.x,p2.y,0.0);
		glEnd();
		glBegin(GL_POINTS);
		glVertex3f(p1.x,p1.y,0.0);
		glVertex3f(p2.x,p2.y,0.0);
		glEnd();
	}

	glDisable(GL_LINE_STIPPLE);

	glPopMatrix();
}

void draw_text_information(bool enable_rot = false)
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0,800.0,0.0,800.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	if(banner_warning_type[cur_banner_id]==0)
		glColor3f(1.0,0.0,0.0);
	else if(banner_warning_type[cur_banner_id]==1)
		glColor3f(0.96,0.77,0.17);
	else
		glColor3f(0.0,1.0,0.0);

	//The main information display
	glRasterPos2i(400,770);
	std::stringstream dban_ss;
	dban_ss << banner_array[cur_banner_id];
	std::string disp_banner = dban_ss.str();
	#ifdef __APPLE__
		_glutBitmapString(GLUT_BITMAP_HELVETICA_18,reinterpret_cast<char const*>(disp_banner.c_str()));
	#else
		glutBitmapString(GLUT_BITMAP_HELVETICA_18,reinterpret_cast<unsigned char const*>(disp_banner.c_str()));
	#endif
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}

void draw_gridlines(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(0.5);

	glColor3f(0.0,0.0,0.0);

	glBegin(GL_LINES);

	for(int i=1;i<21;i++)
	{
		glVertex3f((float)i,0.0,0.0);
		glVertex3f((float)i,20.0,0.0);
		glVertex3f(0.0,(float)i,0.0);
		glVertex3f(20.0,(float)i,0.0);
	}

	glEnd();

	glPopMatrix();
}


void draw_control_points(bool enable_rot = false,bool render_in_3D = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glPointSize( 6.0 );

	glBegin(GL_POINTS);

	if(!render_in_3D)
	{
		if(twodim_axis==0)
		{
			for(int p=0;p<control_point_data_xy.size();++p)
			{
				if(shift_mode_activated && p==clicked_pnt_idx)
					glColor3f(0.0,0.2,0.3);
				else
					glColor3f(0.95,0.2,0.3);
				glVertex3f(control_point_data_xy[p].x,control_point_data_xy[p].y,0.0);
			}
		}
		else
		{
			for(int p=0;p<control_point_data_yz.size();++p)
			{
				if(shift_mode_activated && p==clicked_pnt_idx)
					glColor3f(0.0,0.2,0.3);
				else
					glColor3f(0.95,0.2,0.3);
				glVertex3f(control_point_data_yz[p].x,control_point_data_yz[p].y,0.0);
			}
		}
	}
	else
	{
		for(int p=0;p<control_point_data_xy.size();++p)
		{
			glColor3f(0.95,0.2,0.3);
			glVertex3f(control_point_data_xy[p].x,control_point_data_xy[p].y,0.0);
		}
		for(int p=0;p<control_point_data_yz.size();++p)
		{
			glColor3f(0.95,0.2,0.3);
			glVertex3f(0.0,control_point_data_yz[p].y,control_point_data_yz[p].x);
		}

	}

	glEnd();

	glPopMatrix();
}


void draw_old_control_points(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glPointSize( 6.0 );
	glEnable( GL_POINT_SMOOTH );

	glBegin(GL_POINTS);
	glColor3f(0.0,0.5,0.5);

	for(int p=0;p<repos_points_old.size();++p)
	{
		glVertex3f(repos_points_old[p].x,repos_points_old[p].y,0.0);
	}

	glEnd();
	glDisable( GL_POINT_SMOOTH );
	glPopMatrix();
}

void draw_saved_control_points(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glPointSize( 6.0 );
	glEnable( GL_POINT_SMOOTH );

	glBegin(GL_POINTS);
	glColor3f(0.0,0.5,0.5);

	for(int p=0;p<saved_control_points.size();++p)
	{
		glVertex3f(saved_control_points[p].x,saved_control_points[p].y,0.0);
	}

	glEnd();
	glDisable( GL_POINT_SMOOTH );
	glPopMatrix();
}

void draw_bezier_curve(bool enable_rot = false,bool on_yz_axis=false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<num_points_bezier_curve;++p)
	{
		if(on_yz_axis)
			glVertex3f(0.0,bezier_curve_data[p].y,bezier_curve_data[p].x);
		else
			glVertex3f(bezier_curve_data[p].x,bezier_curve_data[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_old_bezier_curve(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glLineStipple(1,0x3F07);
	glEnable(GL_LINE_STIPPLE);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<num_points_old_bezier_curve;++p)
	{
		glVertex3f(old_bezier_curve_data[p].x,old_bezier_curve_data[p].y,0.0);
	}

	glEnd();

	glDisable(GL_LINE_STIPPLE);

	glPopMatrix();
}

void draw_rational_bezier_curve(bool enable_rot = false,bool on_yz_axis=false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}
	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<num_points_rational_bezier_curve;++p)
	{
		if(on_yz_axis)
			glVertex3f(0.0,rational_bezier_curve_data[p].x,rational_bezier_curve_data[p].y);
		else
			glVertex3f(rational_bezier_curve_data[p].x,rational_bezier_curve_data[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_cubic_bspline_curve(bool enable_rot = false,bool on_yz_axis=false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<num_points_cubicbspline;++p)
	{
		if(on_yz_axis)
			glVertex3f(0.0,cubic_bspline_data[p].y,cubic_bspline_data[p].x);
		else
			glVertex3f(cubic_bspline_data[p].x,cubic_bspline_data[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_subdivision_curves(bool enable_rot = false,bool on_yz_axis=false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<subdivision_curve_data.size();++p)
	{
		if(on_yz_axis)
			glVertex3f(0.0,subdivision_curve_data[p].y,subdivision_curve_data[p].x);
		else
			glVertex3f(subdivision_curve_data[p].x,subdivision_curve_data[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_nurbs_curve(bool enable_rot = false,bool on_yz_axis=false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<nurbs_data.size();++p)
	{
		if(on_yz_axis)
			glVertex3f(0.0,nurbs_data[p].y,nurbs_data[p].x);
		else
			glVertex3f(nurbs_data[p].x,nurbs_data[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_arbbsp_curve(bool enable_rot = false,bool on_yz_axis=false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<arbbsp_data.size();++p)
	{
		if(on_yz_axis)
			glVertex3f(0.0,arbbsp_data[p].y,arbbsp_data[p].x);
		else
			glVertex3f(arbbsp_data[p].x,arbbsp_data[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_sor()
{
	glPushMatrix();


	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );


	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	int counter=0;
	for(int t=0;t<num_quads_sor*4;++t)
	{
		glVertex3f(sor_point_data[sor_indices[t]].x,sor_point_data[sor_indices[t]].y,sor_point_data[sor_indices[t]].z);
	}

	glEnd();

	glPopMatrix();
}

void draw_soe()
{
	glPushMatrix();


	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );


	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	int counter=0;
	for(int t=0;t<num_quads_soe*4;++t)
	{
		glVertex3f(soe_point_data[soe_indices[t]].x,soe_point_data[soe_indices[t]].y,soe_point_data[soe_indices[t]].z);
	}

	glEnd();

	glPopMatrix();
}

void draw_sweep()
{
	glPushMatrix();


	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );


	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	int counter=0;
	for(int t=0;t<num_quads_sweep*4;++t)
	{
		glVertex3f(sweep_point_data[sweep_indices[t]].x,sweep_point_data[sweep_indices[t]].y,sweep_point_data[sweep_indices[t]].z);
	}

	glEnd();

	glPopMatrix();
}

void draw_lofts()
{
	glPushMatrix();


	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );


	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	int counter=0;
	for(int t=0;t<num_quads_loft*4;++t)
	{
		glVertex3f(loft_point_data[loft_indices[t]].x,loft_point_data[loft_indices[t]].y,loft_point_data[loft_indices[t]].z);
	}

	glEnd();

	glPopMatrix();
}

void draw_bezier_surface()
{
	glPushMatrix();

	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	for(int i=0;i<num_samples_in_u-1;i++)
	{
		for(int j=0;j<num_samples_in_v-1;j++)
		{
			glVertex3f(bezier_surface_samples[i][j].x,bezier_surface_samples[i][j].y,bezier_surface_samples[i][j].z);
			glVertex3f(bezier_surface_samples[i+1][j].x,bezier_surface_samples[i+1][j].y,bezier_surface_samples[i+1][j].z);
			glVertex3f(bezier_surface_samples[i+1][j+1].x,bezier_surface_samples[i+1][j+1].y,bezier_surface_samples[i+1][j+1].z);
			glVertex3f(bezier_surface_samples[i][j+1].x,bezier_surface_samples[i][j+1].y,bezier_surface_samples[i][j+1].z);
		}
	}

	glEnd();

	glPopMatrix();
}

void draw_bspline_surface()
{
	glPushMatrix();

	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	for(int i=0;i<num_samples_in_u_bsp-1;i++)
	{
		for(int j=0;j<num_samples_in_v_bsp-1;j++)
		{
			glVertex3f(bspline_surface_samples[i][j].x,bspline_surface_samples[i][j].y,bspline_surface_samples[i][j].z);
			glVertex3f(bspline_surface_samples[i+1][j].x,bspline_surface_samples[i+1][j].y,bspline_surface_samples[i+1][j].z);
			glVertex3f(bspline_surface_samples[i+1][j+1].x,bspline_surface_samples[i+1][j+1].y,bspline_surface_samples[i+1][j+1].z);
			glVertex3f(bspline_surface_samples[i][j+1].x,bspline_surface_samples[i][j+1].y,bspline_surface_samples[i][j+1].z);
		}
	}

	glEnd();

	glPopMatrix();
}

void draw_trimmednurbs_surface()
{
	glPushMatrix();

	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	glBegin(GL_QUADS);

	glColor3f(0.0,0.0,0.0);

	for(int i=0;i<tnurbs_samples_in_x-1;i++)
	{
		for(int j=0;j<tnurbs_samples_in_y-1;j++)
		{
			glVertex3f(tnurbs_surface_samples[i][j].x,tnurbs_surface_samples[i][j].y,tnurbs_surface_samples[i][j].z);
			glVertex3f(tnurbs_surface_samples[i+1][j].x,tnurbs_surface_samples[i+1][j].y,tnurbs_surface_samples[i+1][j].z);
			glVertex3f(tnurbs_surface_samples[i+1][j+1].x,tnurbs_surface_samples[i+1][j+1].y,tnurbs_surface_samples[i+1][j+1].z);
			glVertex3f(tnurbs_surface_samples[i][j+1].x,tnurbs_surface_samples[i][j+1].y,tnurbs_surface_samples[i][j+1].z);
		}
	}

	glEnd();

	glPopMatrix();
}

void draw_history_curve_on_xy()
{
	glPushMatrix();

	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(1.0);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<history_curve.size();++p)
	{
		glVertex3f(history_curve[p].x,history_curve[p].y,0.0);
	}

	glEnd();

	glPopMatrix();
}

void draw_control_polyhedron()
{
	glPushMatrix();


	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(1.0);

	glLineStipple(1,0x3F07);
	glEnable(GL_LINE_STIPPLE);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	glBegin(GL_QUADS);

	glColor3f(1.0,0.0,0.0);

	for(int i=0;i<m-1;i++)
	{
		for(int j=0;j<n-1;j++)
		{
			glVertex3f(control_net_array[i][j].x,control_net_array[i][j].y,control_net_array[i][j].z);
			glVertex3f(control_net_array[i+1][j].x,control_net_array[i+1][j].y,control_net_array[i+1][j].z);
			glVertex3f(control_net_array[i+1][j+1].x,control_net_array[i+1][j+1].y,control_net_array[i+1][j+1].z);
			glVertex3f(control_net_array[i][j+1].x,control_net_array[i][j+1].y,control_net_array[i][j+1].z);
		}
	}

	glEnd();

	glDisable(GL_LINE_STIPPLE);

	glPopMatrix();
}

void draw_control_points_in_3d()
{
	glPushMatrix();

	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glPointSize( 6.0 );

	glColor3f(0.95,0.2,0.3);

	glBegin(GL_POINTS);

	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			glVertex3f(control_net_array[i][j].x,control_net_array[i][j].y,control_net_array[i][j].z);
		}
	}
	glEnd();

	glPopMatrix();
}

void draw_saved_curve_on_xy(bool enable_rot = false)
{
	glPushMatrix();

	if(enable_rot)
	{
		glRotatef(x_angle,0,1,0);
		glRotatef(y_angle,1,0,0);
		glScalef(scale_size,scale_size,scale_size);
	}

	glLineWidth(1.0);

	glLineStipple(1,0x3F07);
	glEnable(GL_LINE_STIPPLE);

	glBegin(GL_LINE_STRIP);

	glColor3f(0.0,0.0,0.0);

	for(int p=0;p<saved_curve_samples.size();++p)
	{
		glVertex3f(saved_curve_samples[p].x,saved_curve_samples[p].y,0.0);
	}

	glEnd();

	glDisable(GL_LINE_STIPPLE);

	glPopMatrix();
}

void draw_mesh(Mesh* m,std::vector<glm::vec3> &colvec)
{
	glPushMatrix();


	glRotatef(x_angle,0,1,0);
	glRotatef(y_angle,1,0,0);
	glScalef(scale_size,scale_size,scale_size);

	glLineWidth(0.5);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glBegin(GL_TRIANGLES);

	glColor3f(1.0,0.0,0.0);

	int n_faces = m->GetNumberFacets();

	for(int i=0;i<n_faces;i++)
	{
		glm::vec3 facecolor = colvec[i];
		glColor3f(facecolor.x,facecolor.y,facecolor.z);
		// TopoFacet ithface = m->GetFacet(face_render_cntr);
		TopoFacet ithface = m->GetFacet(i);
		int facet_degree = ithface.GetNumberVertices();
		// std::cout << facet_degree << "\n";
		for(int t=0;t<facet_degree-2;t++)
		{
			int i1 = ithface.GetVertexInd(0);
			int i2 = ithface.GetVertexInd(t+1);
			int i3 = ithface.GetVertexInd(t+2);
			GeomVert v1 = m->GetGeomVertex(i1);
			glVertex3f(v1.GetCo(0),v1.GetCo(1),v1.GetCo(2));
			GeomVert v2 = m->GetGeomVertex(i2);
			glVertex3f(v2.GetCo(0),v2.GetCo(1),v2.GetCo(2));
			GeomVert v3 = m->GetGeomVertex(i3);
			glVertex3f(v3.GetCo(0),v3.GetCo(1),v3.GetCo(2));
		}
	}


	glEnd();

	glPopMatrix();
}

void display()
{
	//Split into 2 portions
	//This one is 2D
	glEnable(GL_DEPTH_TEST);
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glViewport(0,0,800,800);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60,1.0,0.1,100);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(10.5,10.6,20,10.5,10.6,0,0,1,0);
	draw_axes(false,twodim_axis);
	if(enable_grid_lines)
	{
		draw_gridlines();
	}
	if(enable_control_polygon)
	{
		draw_control_polygon();
	}
	draw_control_points();
	draw_residues();
	switch(curve_type)
	{
		case 0 : draw_bezier_curve(); break;
		case 1 : draw_cubic_bspline_curve(); break;
		case 2 : draw_subdivision_curves(); break;
		case 3 : draw_rational_bezier_curve(); break;
		case 4 : if(arbbsp_ready_for_rendering) draw_arbbsp_curve(); break;
		case 5 : if(nurbs_ready_for_rendering) draw_nurbs_curve(); break;
	}
	if(add_disp_c1_cont)
	{
		draw_old_control_points();
		draw_old_bezier_curve();
	}
	if(lofting_context_mode)
	{
		draw_saved_control_points(false);
		draw_saved_curve_on_xy();
	}
	draw_text_information();

	//Now the 3D view
	glViewport(801,0,800,800);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60,1.0,0.1,100);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(30,30,30,0,0,0,0,1,0);
	draw_axes_3d(true);
	if(ds_exists)
	{
		draw_mesh(subdiv_res_ds,ds_facet_colors);
	}
	else if(cc_exists)
	{
		draw_mesh(subdiv_res_cc,cc_facet_colors);
	}
	else if(loop_exists)
	{
		draw_mesh(subdiv_res_loop,loop_facet_colors);
	}
	else if(off_exists)
	{
		draw_mesh(shape,shape_facet_colors);
	}
	else if(bez_surface_rendered)
	{
		draw_bezier_surface();
		if(enable_control_polyhedron)
		{
			draw_control_points_in_3d();
			draw_control_polyhedron();
		}
	}
	else if(bsp_surface_rendered)
	{
		draw_bspline_surface();
		if(enable_control_polyhedron)
		{
			draw_control_points_in_3d();
			draw_control_polyhedron();
		}
	}
	else if(tnurbs_surface_rendered)
	{
		draw_trimmednurbs_surface();
		if(enable_control_polyhedron)
		{
			draw_control_points_in_3d();
			draw_control_polyhedron();
		}
	}
	else if(sor_display)
	{
		draw_sor();
	}
	else if(soe_display)
	{
		draw_soe();
	}
	else if(sweep_display)
	{
		//Actual Curve
		draw_sweep();
		draw_control_points(true,true);
		draw_history_curve_on_xy(); //The generator wire
		switch(curve_type)
		{
			case 0 : draw_bezier_curve(true,true); break;
			case 1 : draw_cubic_bspline_curve(true,true); break;
			case 2 : draw_subdivision_curves(true,true); break;
			case 3 : draw_rational_bezier_curve(true,true); break;
			case 4 : if(arbbsp_ready_for_rendering) draw_arbbsp_curve(true,true); break;
			case 5 : if(nurbs_ready_for_rendering) draw_nurbs_curve(true,true); break;
		}

	}
	else if(loft_display)
	{
		draw_lofts();
	}
	glutSwapBuffers();
}

glm::vec2 get_grid_position(int x, int y)
{
	glm::vec2 res;

	float xf = (float)x;
	float yf = (float)y;

	float x_offset = (xf-37.0);
	float x_gap = (730.0-37.0)/20.0;
	float x_coods = x_offset/x_gap;

	float y_offset = (770.0-yf);
	float y_gap = (770.0-76.0)/20.0;
	float y_coods = y_offset/y_gap;

	res.x = x_coods;
	res.y = y_coods;

	return res;
	// std::cout << x_coods << "|" << y_coods << "\n";
}

void mouse_func(int btn, int state, int x, int y)
{
	glm::vec2 pnt_vec;
    if(state==GLUT_DOWN)
    {
    	if(btn==GLUT_LEFT_BUTTON)
    	{
    		press_x = x;
    		press_y = y;
    		clicked_type = CLICK_MODE;
    		if(x>=37 && x<=730 && y>=76 && y<=770 )
	    	{
	    		add_disp_c1_cont = false;
	    		pnt_vec = get_grid_position(x,y);
	    		if(twodim_axis==0)
	    			control_point_data_xy.push_back(pnt_vec);
	    		else
	    			control_point_data_yz.push_back(pnt_vec);
	    		Point_3 cgal_p(pnt_vec.x,pnt_vec.y,0.0);
	    		if(io_type==1)
	    		{
	    			if(control_point_data_xy.size()>2)
	    				if(recons_type==0)
	    					curve_reconstruction_crust();
	    				else
								curve_reconstruction_NN_crust();
	    		}
				switch(curve_type)
    			{
    				case 0 : generate_bezier_curve(); break;
    				case 1 : generate_cubicbspline_curve(); break;
    				case 2 : generate_subdivision_curves(); break;
    				case 3 : generate_rational_bezier_curve(); break;
    			}
	    	}
	    	else if(x>=801)
	    	{
	    		xform_mode = XFORM_ROTATE;
	    	}
	    	else
	    	{
	    		cur_banner_id = 3;
	    		std::cout << "[WARNING] Mouse clicked outside grid. Will not be registered.\n";
	    	}
	    }
	    else if(btn==GLUT_RIGHT_BUTTON)
	    {
	    	if(x>=37 && x<=730 && y>=76 && y<=770 )
	    	{
	    		pnt_vec = get_grid_position(x,y);
	    		if(twodim_axis==0)
	    		{
		    		for(int i=0;i<control_point_data_xy.size();++i)
		    		{
		    			float dist = glm::length(pnt_vec-control_point_data_xy[i]);
		    			if(dist<0.2)
		    			{
		  					clicked_pnt_idx = i;
		  					shift_mode_activated = true;
		  					clicked_type = SHIFT_MODE;
		    				break;
		    			}
		    		}
	    		}
	    		else
	    		{
	    			for(int i=0;i<control_point_data_yz.size();++i)
		    		{
		    			float dist = glm::length(pnt_vec-control_point_data_yz[i]);
		    			if(dist<0.2)
		    			{
		  					clicked_pnt_idx = i;
		  					shift_mode_activated = true;
		  					clicked_type = SHIFT_MODE;
		    				break;
		    			}
		    		}

	    		}
	    	}
	    	else if(x>=801)
	    	{
	    		xform_mode = XFORM_SCALE;
	    	}
	    	else
	    	{
	    		cur_banner_id = 3;
	    		std::cout << "[WARNING] Mouse clicked outside grid. Will not be registered.\n";
	    	}
	    }
    }
    else if(state==GLUT_UP)
    {
    	xform_mode = XFORM_NONE;
    	if(btn==GLUT_RIGHT_BUTTON)
    	{
    		if(shift_mode_activated)
    		{
    			if(io_type==1)
	    		{
	    			if(control_point_data_xy.size()>2)
	    				if(recons_type==0)
	    					curve_reconstruction_crust();
	    				else
								curve_reconstruction_NN_crust();
	    		}
    			switch(curve_type)
	    		{
	    			case 0 : generate_bezier_curve(true); break;
	    			case 1 : generate_cubicbspline_curve(true); break;
	    			case 2 : generate_subdivision_curves(); break;
	    			case 3 : generate_rational_bezier_curve(); break;
	    		}
    		}
    		clicked_type = NONE;
    		clicked_pnt_idx = -1;
    		shift_mode_activated = false;
    	}
    }
    add_disp_c1_cont = false;
    glutPostRedisplay();
}

void process_duplicity()
{
	int x,y;
	x = global_mouse_x;
	y = global_mouse_y;
	// glfwGetMousePos(&x,&y);
	glm::vec2 pnt_vec;
	bool point_clicked = false;
	if(x>=37 && x<=730 && y>=76 && y<=770 )
	{
		pnt_vec = get_grid_position(x,y);
		for(int i=0;i<control_point_data_xy.size();++i)
		{
			float dist = glm::length(pnt_vec-control_point_data_xy[i]);
			if(dist<0.2)
			{
				glm::vec2 data = control_point_data_xy[i];
				control_point_data_xy.insert(control_point_data_xy.begin()+i,data);
				std::cout << i << " detected for duplicity\n";
				//Insert an add fot duplicity
				//Redraw the curve
				switch(curve_type)
				{
					case 0 : generate_bezier_curve(); break;
					case 1 : generate_cubicbspline_curve(); break;
					case 2 : generate_subdivision_curves(); break;
					case 3 : generate_rational_bezier_curve(); break;
				}
				cur_banner_id = 11;
				point_clicked = true;
				break;
			}
		}
		if(!point_clicked)
			cur_banner_id = 12;
	}
	else
	{
		cur_banner_id = 3;
	}
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void process_point_deletion()
{
	int x,y;
	x = global_mouse_x;
	y = global_mouse_y;
	// glfwGetMousePos(&x,&y);
	glm::vec2 pnt_vec;
	bool point_clicked = false;
	if(x>=37 && x<=730 && y>=76 && y<=770 )
	{
		pnt_vec = get_grid_position(x,y);
		for(int i=0;i<control_point_data_xy.size();++i)
		{
			float dist = glm::length(pnt_vec-control_point_data_xy[i]);
			if(dist<0.2)
			{
				control_point_data_xy.erase(control_point_data_xy.begin()+i);
				std::cout << i << " detected for duplicity\n";
				//Insert an add fot duplicity
				//Redraw the curve
				if(io_type==1)
				{
					if(control_point_data_xy.size()>2)
						if(recons_type==0)
							curve_reconstruction_crust();
						else
							curve_reconstruction_NN_crust();
				}
				switch(curve_type)
				{
					case 0 : generate_bezier_curve(); break;
					case 1 : generate_cubicbspline_curve(); break;
					case 2 : generate_subdivision_curves(); break;
					case 3 : generate_rational_bezier_curve(); break;
				}
				cur_banner_id = 35;
				point_clicked = true;
				break;
			}
		}
		if(!point_clicked)
			cur_banner_id = 12;
	}
	else
	{
		cur_banner_id = 3;
	}
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void analyze_mesh(Mesh* m)
{
	std::cout << "Number of vertices : " << m->GetNumberVertices() << "\n";
	std::cout << "NUmber of edges : " << m->GetNumberEdges() << "\n";
	std::cout << "Number of faces : " << m->GetNumberFacets() << "\n";
}

void export_mesh_to_obj(Mesh* obj,char* f_name)
{
	std::ofstream file_output(f_name);
	int n_vertices = obj->GetNumberVertices();
	for(int i=0;i<n_vertices;i++)
	{
		GeomVert gv = obj->GetGeomVertex(i);
		file_output << "v " << gv.GetCo(0) << " " << gv.GetCo(1) << " " << gv.GetCo(2) << "\n";
	}
	int n_faces = obj->GetNumberFacets();
	for(int i=0;i<n_faces;i++)
	{
		TopoFacet f = obj->GetFacet(i);
		int n_sur_vertices = f.GetNumberVertices();
		file_output << "f " ;
		for(int k=0;k<n_sur_vertices;k++)
		{
			file_output << (f.GetVertexInd(k)+1) << " ";
		}
		file_output << "\n";
	}
	file_output.close();
}


void keybd_func(unsigned char key,int x,int y)
{
	switch(key)
	{
		case 'd' : case 'D' :
			process_duplicity(); break;
		case 'x' : case 'X' :
			process_point_deletion(); break;
		case 'q' : case 'Q' :
			exit(0);
		case 's' : case 'S' :
			if(off_exists)
				analyze_mesh(shape);
			break;
		case 'z' :
			face_render_cntr++;
			std::cout << "Value of debug counter is" << face_render_cntr << "\n";
			break;
	}

	glutPostRedisplay();
}

void motion_func(int x, int y)
{
	if (xform_mode==XFORM_ROTATE) {

      x_angle += (x - press_x)/5.0;
      if (x_angle > 180) x_angle -= 360;
      else if (x_angle <-180) x_angle += 360;
      press_x = x;

      y_angle += (y - press_y)/5.0;
      if (y_angle > 180) y_angle -= 360;
      else if (y_angle <-180) y_angle += 360;
      press_y = y;
     }
	else if (xform_mode == XFORM_SCALE){
      float old_size = scale_size;
      scale_size *= (1+ (y - press_y)/60.0);
      if (scale_size <0) scale_size = old_size;
      press_y = y;
    }
	if(clicked_type==SHIFT_MODE)
	{
		if(clicked_pnt_idx !=-1)
		{
			if(twodim_axis==0)
				control_point_data_xy[clicked_pnt_idx] = get_grid_position(x,y);
			else
				control_point_data_yz[clicked_pnt_idx] = get_grid_position(x,y);
		}

		switch(curve_type)
		{
			case 0 : generate_bezier_curve(true); break;
			case 1 : generate_cubicbspline_curve(true); break;
			case 2 : generate_subdivision_curves(); break;
			case 3 : generate_rational_bezier_curve(); break;
		}
	}
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void passive_motion_func(int x,int y)
{
	global_mouse_x = x;
	global_mouse_y = y;
}

void simple_redisplay_cb_func(int id)
{
	glutPostRedisplay();
}

void bezier_close_cb(int id)
{
	if(curve_type!=0)
		return;

	glm::vec2 p0,pnm1;

	if(twodim_axis==0)
	{
		p0 = control_point_data_xy.front();
		pnm1 = control_point_data_xy.back();
	}
	else
	{
		p0 = control_point_data_yz.front();
		pnm1 = control_point_data_yz.back();
	}

	//Make them collinear now
	glm::vec2 v = pnm1 - p0;
	glm::vec2 p1 = p0 - v;

	if(twodim_axis==0)
	{
		control_point_data_xy[1] = p1;
		control_point_data_xy.push_back(p0);
	}
	else
	{
		control_point_data_yz[1] = p1;
		control_point_data_yz.push_back(p0);
	}
	possibly_closed = true;
	generate_bezier_curve();
	glutPostRedisplay();
}

void spinner_cb_func(int id)
{
	//Add the refinement code here
	if(curve_type==2)
	{
		generate_subdivision_curves();
		add_disp_c1_cont = false;
		glutPostRedisplay();
	}
}

void reset_cb_func(int id)
{
	//Clear control point vector
	control_point_data_xy.clear();
	control_point_data_yz.clear();
	repos_points_old.clear();

	//set all the rendered points value to 0
	num_points_old_bezier_curve = 0;
	num_points_bezier_curve=0;
	num_points_rational_bezier_curve=0;
	num_points_cubicbspline=0;

	//Set the subdivision curve deque to 0
	subdivision_curve_data.clear();

	//Set the rational weight value to 0
	number_of_weights = 0;
	number_of_weights_nurbs = 0;

	//Banner id reset to 1
	cur_banner_id = 0;

	nurbs_data.clear();
	nurbs_ready_for_rendering = false;

	arbbsp_data.clear();
	arbbsp_ready_for_rendering = false;

	add_disp_c1_cont = false;

	//Clear sor data if any
	sor_point_data.clear();
	sor_indices.clear();
	num_quads_sor=0;
	sor_display  = false;

	//Clear soe data if any
	soe_point_data.clear();
	soe_indices.clear();
	num_quads_soe=0;
	soe_display  = false;

	history_curve.clear();
	sweep_point_data.clear();
	sweep_indices.clear();
	num_quads_sweep=0;
	sweep_display  = false;

	lofting_context_mode = false;
	saved_curve_samples.clear();
	saved_control_points.clear();
	loft_display = false;
	num_quads_loft = 0;
	loft_point_data.clear();
	loft_indices.clear();

	if(off_exists)
	{
		delete shape;
		off_exists = false;
		shape_facet_colors.clear();
	}
	if(ds_exists)
	{
		delete subdiv_res_ds;
		ds_exists = false;
		ds_facet_colors.clear();
	}
	if(cc_exists)
	{
		delete subdiv_res_cc;
		cc_exists = false;
		cc_facet_colors.clear();
	}
	if(loop_exists)
	{
		delete subdiv_res_loop;
		loop_exists = false;
		loop_facet_colors.clear();
	}

	m=n=0;
	possibly_closed = false;
	bez_surface_rendered = false;
	num_samples_in_u_bsp = 0;
	num_samples_in_v_bsp = 0;
	bsp_surface_rendered = false;

	glutPostRedisplay();
}

void samples_cb_func(int id)
{
	if(curve_type==0)
	{
		generate_bezier_curve();
		add_disp_c1_cont = false;
		glutPostRedisplay();
	}
}

void radio_button_cb_func(int id)
{
	switch(curve_type)
	{
		case 0 : generate_bezier_curve(); break;
		case 1 : generate_cubicbspline_curve(); break;
		case 2 : generate_subdivision_curves(); break;
		case 3 : generate_rational_bezier_curve(); break;
		case 4 : if(arbbsp_ready_for_rendering) cur_banner_id = 20; break;
		case 5 : if(nurbs_ready_for_rendering) cur_banner_id = 17; break;
	}
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void twodim_axis_radio_cb_func(int id)
{
	//Save the old points in the history data structure
	history_curve.clear();
	if(curve_type==0)
	{
		for(int i=0;i<num_subdivisions+1;i++)
		{
			history_curve.push_back(bezier_curve_data[i]);
		}
		num_points_bezier_curve = 0;

	}
	else if(curve_type==1)
	{
		for(int i=0;i<num_points_cubicbspline;i++)
		{
			history_curve.push_back(cubic_bspline_data[i]);
		}
		num_points_cubicbspline = 0;
	}
	else if(curve_type==2)
	{
		for(int i=0;i<subdivision_curve_data.size();i++)
		{
			history_curve.push_back(subdivision_curve_data[i]);
		}
		subdivision_curve_data.clear();
	}
	else if(curve_type==3)
	{
		for(int i=0;i<num_points_rational_bezier_curve;i++)
		{
			history_curve.push_back(rational_bezier_curve_data[i]);
		}
		num_points_rational_bezier_curve = 0;
	}
	else if(curve_type==4)
	{
		history_curve = arbbsp_data;
		arbbsp_data.clear();
	}
	else if(curve_type==5)
	{
		history_curve = nurbs_data;
		nurbs_data.clear();
	}
	switch(curve_type)
	{
		case 0 : generate_bezier_curve(); break;
		case 1 : generate_cubicbspline_curve(); break;
		case 2 : generate_subdivision_curves(); break;
		case 3 : generate_rational_bezier_curve(); break;
		case 4 : if(arbbsp_ready_for_rendering) cur_banner_id = 20; break;
		case 5 : if(nurbs_ready_for_rendering) cur_banner_id = 17; break;
	}
	if(twodim_axis==0)
		cur_banner_id = 22;
	else
		cur_banner_id = 21;
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void load_rbsweight_cb_func(int id)
{
	std::ifstream in_file("rational_bezier_weights.dat");
	std::string line;
	std::string weight_str;
	int counter=0;
	while (true)
	{
    	float weight;
    	in_file >> weight;
    	weight_rational_bezier[counter++] = weight;
    	if( in_file.eof() )
    		break;
	}
	number_of_weights = counter;
	cur_banner_id = 8;
	in_file.close();
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void load_knotvectors_cb_func(int id)
{
	std::ifstream in_file("knot_vector.dat");
	std::string line;
	int counter=0;
	while (true)
	{
    	float kv;
    	in_file >> kv;
    	arbbsp_knot_vector[counter++] = kv;
    	if( in_file.eof() )
    		break;
	}
	number_of_knot_vectors = counter;
	cur_banner_id = 18;
	in_file.close();
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

/*
void load_knotvectorsx2_cb_func(char* file_name,int id)
{
	std::ifstream in_file(file_name);
	std::string line;
	int counter=0;
	while (true)
	{
    	float kv;
    	in_file >> kv;
    	if(id==0)
    		knot_vector_i[counter++] = kv;
    	else
    		knot_vector_j[counter++] = kv;
    	if( in_file.eof() )
    		break;
	}
	if(id==0)
		number_of_knot_vectors_i = counter;
	else
		number_of_knot_vectors_j = counter;
}
*/

void load_nurbweight_cb_func(int id)
{
	std::ifstream in_file("nurbs_weights.dat");
	std::string line;
	std::string weight_str;
	int counter=0;
	while (true)
	{
    	float weight;
    	in_file >> weight;
    	weight_nurbs[counter++] = weight;
    	if( in_file.eof() )
    		break;
	}
	number_of_weights_nurbs = counter;
	cur_banner_id = 15;
	in_file.close();
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void c1_split_cb(int id)
{
	int n_org = control_point_data_xy.size();
	//check if its of the form 3k+1
	//if not remove it to make it of the 3k+1 form
	int prob_k = (n_org-1)/3;
	if(prob_k<2)
	{
		//Show error and return
		cur_banner_id = 14;
		glutPostRedisplay();
		return;
	}
	//Draw the old curve
	int counter = 0;
	num_points_old_bezier_curve=0;
	for(int i=0;i<prob_k;i++)
	{
		std::vector<glm::vec2> res_k_seg = get_points_on_cubic_bezier(control_point_data_xy[counter],control_point_data_xy[counter+1],control_point_data_xy[counter+2],control_point_data_xy[counter+3]);
		for(int j=0;j<11;j++)
		{
			old_bezier_curve_data[num_points_old_bezier_curve++] = res_k_seg[j];
		}
		counter +=3;
	}
	//Reposition CP
	for(int i=0;i<prob_k-1;i++)
	{
		int idx_to_mod = 3*(i+1);
		repos_points_old.push_back(control_point_data_xy[idx_to_mod+1]);
		glm::vec2 vecp3p2 = control_point_data_xy[idx_to_mod]-control_point_data_xy[idx_to_mod-1];
		control_point_data_xy[idx_to_mod+1] = control_point_data_xy[idx_to_mod] + vecp3p2;
	}
	//Generate new bezier curve
	counter = 0;
	num_points_bezier_curve=0;
	for(int i=0;i<prob_k;i++)
	{
		std::vector<glm::vec2> res_k_seg = get_points_on_cubic_bezier(control_point_data_xy[counter],control_point_data_xy[counter+1],control_point_data_xy[counter+2],control_point_data_xy[counter+3]);
		for(int j=0;j<11;j++)
		{
			bezier_curve_data[num_points_bezier_curve++] = res_k_seg[j];
		}
		counter +=3;
	}
	//Render
	cur_banner_id = 13;
	add_disp_c1_cont = true;
	glutPostRedisplay();
}

void gen_nurbs_cb(int id)
{
	int num_cp = control_point_data_xy.size();
	if(number_of_weights_nurbs < num_cp )
	{
		cur_banner_id = 16;
		glutPostRedisplay();
		return;
	}
	int n = num_cp-1;
	float* knot_vec = generate_nonuniform_knot_vector(nurb_k,n);
	int max_val_u = n-nurb_k+2;
	// std::cout << max_val_u << "\n";
	float u = 0.0;
	float gap = 0.05;
	nurbs_data.clear();
	while(u<1.05)// while(u<=(float)(max_val_u))|| abs(u-1.0)<TOLERANCE)
	{
		// std::cout << "Processing u-val:" << u << "\n";
		float* c_vec = generate_coefficients_for_u(nurb_k,n,knot_vec,u);
		glm::vec2 output(0.0,0.0);
		float denom = 0.0;
		for(int i=0;i<=n;i++)
		{
			output += (weight_nurbs[i]*c_vec[i]*control_point_data_xy[i]);
			denom += (weight_nurbs[i]*c_vec[i]);
		}
		float denom_inv = (float)1.0/denom;
		output = denom_inv * output ;
		nurbs_data.push_back(output);
		u +=gap;
	}
	cur_banner_id = 17;
	nurbs_ready_for_rendering = true;
	glutPostRedisplay();
}

void gen_arbbsp_cb(int id)
{
	//First do a sanity check
	int n_org = control_point_data_xy.size();
	int n = n_org-1;
	if(number_of_knot_vectors<(n+arbbsp_k+1))
	{
		cur_banner_id = 19;
		add_disp_c1_cont = false;
		glutPostRedisplay();
	}
	//Sort the array
	std::sort(arbbsp_knot_vector,arbbsp_knot_vector+number_of_knot_vectors,std::greater<float>());
	//Find the max and min
	float kv_min = arbbsp_knot_vector[0];
	float kv_max = arbbsp_knot_vector[number_of_knot_vectors-1];
	float norm_denom = kv_max-kv_min;
	for(int i=0;i<number_of_knot_vectors;i++)
	{
		//@TODO: Check for negative numbers or something
		arbbsp_knot_vector[i] = (arbbsp_knot_vector[i]-kv_min)/norm_denom;
		std::cout << arbbsp_knot_vector[i] << " ";
	}
	std::cout << "\n";
	arbbsp_data.clear();
	float u =0.0;
	float gap = 0.05;
	while(u<1.05)
	{
		float* c_vec = generate_coefficients_for_u(arbbsp_k,n,arbbsp_knot_vector,u);
		glm::vec2 output(0.0,0.0);
		float denom = 0.0;
		for(int i=0;i<=n;i++)
		{
			output += (c_vec[i]*control_point_data_xy[i]);
		}
		arbbsp_data.push_back(output);
		u +=gap;
	}
	cur_banner_id = 20;
	arbbsp_ready_for_rendering = true;
	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void draw_next_curve(int id)
{
	//First save the current type of rendered curve in the
}

Mesh* generate_mesh_datastructure(std::vector<glm::vec3> &pointset,std::vector<int> &indices,int num_quads)
{
	Mesh* result = new Mesh();

	for(int i=0;i<num_quads*4;i=i+4)
	{
		std::vector<GeomVert> temp_facet;

		glm::vec3 v1 = pointset[indices[i+0]];
		temp_facet.push_back(GeomVert(v1.x,v1.y,v1.z));
		glm::vec3 v2 = pointset[indices[i+1]];
		temp_facet.push_back(GeomVert(v2.x,v2.y,v2.z));
		glm::vec3 v3 = pointset[indices[i+2]];
		temp_facet.push_back(GeomVert(v3.x,v3.y,v3.z));
		glm::vec3 v4 = pointset[indices[i+3]];
		temp_facet.push_back(GeomVert(v4.x,v4.y,v4.z));

		result->AddFacet(temp_facet);
	}

	return result;
}

void bezier_surface_sample_change_cb(int placeholder)
{
	if(bez_surface_rendered)
	{
		generate_bezier_surface();
		glutPostRedisplay();
	}
}

void generate_surface_cb(int id)
{
	if(surface_type==0)
	{
		generate_bezier_surface();
		cur_banner_id = 25;
	}
	else if(surface_type==1)
	{
		generate_bspline_surface();
		cur_banner_id = 26;
	}
	else
	{
		generate_trimmednurbs_surface();
	}
	glutPostRedisplay();
}

void read_off_cb(int id)
{
	if(off_exists)
	{
		delete shape;
		shape = new Mesh();
		off_exists = false;
	}
	else
	{
		shape = new Mesh();
	}
	bool res = process_off_file("Off_Files/octa.off",&shape);
	if(res)
	{
		cur_banner_id = 27;
		off_exists = true;
	}
	else
	{
		cur_banner_id = 28;
		off_exists = false;
		delete shape;
	}
	generate_face_colors(shape,shape_facet_colors);
	glutPostRedisplay();
}

void export_shape_obj_cb(int id)
{
	if(ds_exists)
	{
		export_mesh_to_obj(subdiv_res_ds,"Obj_Files/export_ds.obj");
	}
	else if(cc_exists)
	{
		export_mesh_to_obj(subdiv_res_cc,"Obj_Files/export_cc.obj");
	}
	else if(loop_exists)
	{
		export_mesh_to_obj(subdiv_res_loop,"Obj_Files/export_loop.obj");
	}
}

bool check_mesh_sanity_for_loop(Mesh* input)
{
	int n_faces = input->GetNumberFacets();
	for(int i=0;i<n_faces;i++)
	{
		TopoFacet f = input->GetFacet(i);
		int n_v = f.GetNumberVertices();
		if(n_v!=3)
			return false;
	}
	return true;
}

void subdivision_processing_cb(int id)
{
	if(off_exists)
	{
		if(subdiv_type==0)
		{
			if(!ds_exists)
			{
				subdiv_res_ds = new Mesh();
				doo_sabine(shape,&subdiv_res_ds);
				ds_exists = true;
				generate_face_colors(subdiv_res_ds,ds_facet_colors);
			}
			else
			{
				Mesh* new_mesh = new Mesh();
				doo_sabine(subdiv_res_ds,&new_mesh);
				delete subdiv_res_ds;
				subdiv_res_ds = new_mesh;
				generate_face_colors(subdiv_res_ds,ds_facet_colors);
			}
			cur_banner_id = 30;
		}
		else if(subdiv_type==1)
		{
			if(!cc_exists)
			{
				subdiv_res_cc = new Mesh();
				catmull_clark(shape,&subdiv_res_cc);
				cc_exists = true;
				generate_face_colors(subdiv_res_cc,cc_facet_colors);
			}
			else
			{
				Mesh* new_mesh = new Mesh();
				catmull_clark(subdiv_res_cc,&new_mesh);
				delete subdiv_res_cc;
				subdiv_res_cc = new_mesh;
				generate_face_colors(subdiv_res_cc,cc_facet_colors);
			}
			cur_banner_id = 31;
		}
		else if(subdiv_type==2)
		{
			if(!loop_exists)
			{
				bool ok_to_proceed = check_mesh_sanity_for_loop(shape);
				if(!ok_to_proceed)
				{
					cur_banner_id = 33;
					glutPostRedisplay();
					return;
				}
				subdiv_res_loop = new Mesh();
				loop(shape,&subdiv_res_loop);
				loop_exists = true;
				generate_face_colors(subdiv_res_loop,loop_facet_colors);
			}
			else
			{
				Mesh* new_mesh = new Mesh();
				loop(subdiv_res_loop,&new_mesh);
				delete subdiv_res_loop;
				subdiv_res_loop = new_mesh;
				generate_face_colors(subdiv_res_loop,loop_facet_colors);
			}
			cur_banner_id = 32;
		}
	}
	else
	{
		cur_banner_id = 29;
	}
	glutPostRedisplay();
}

void switch_drawing_context(int id)
{
	lofting_context_mode = true;
	saved_control_points.clear();
	saved_curve_samples.clear();
	saved_curve_type = -1;
	//Saved curve sample
	saved_control_points = control_point_data_xy;
	if(curve_type==0)
	{
		for(int i=0;i<num_subdivisions+1;i++)
		{
			saved_curve_samples.push_back(bezier_curve_data[i]);
		}
		num_points_bezier_curve = 0;
	}
	else if(curve_type==1)
	{
		for(int i=0;i<num_points_cubicbspline;i++)
		{
			saved_curve_samples.push_back(cubic_bspline_data[i]);
		}
		num_points_cubicbspline = 0;
	}
	else if(curve_type==2)
	{
		for(int i=0;i<subdivision_curve_data.size();i++)
		{
			saved_curve_samples.push_back(subdivision_curve_data[i]);
		}
		subdivision_curve_data.clear();
	}
	else if(curve_type==3)
	{
		for(int i=0;i<num_points_rational_bezier_curve;i++)
		{
			saved_curve_samples.push_back(rational_bezier_curve_data[i]);
		}
		num_points_rational_bezier_curve = 0;
	}
	else if(curve_type==4)
	{
		saved_curve_samples = arbbsp_data;
		arbbsp_data.clear();
		arbbsp_ready_for_rendering = false;
	}
	else if(curve_type==5)
	{
		saved_curve_samples = nurbs_data;
		nurbs_data.clear();
		nurbs_ready_for_rendering = false;
	}
	saved_curve_type = curve_type;
	// Erase the sample points
	control_point_data_xy.clear();

	cur_banner_id = 24;

	add_disp_c1_cont = false;
	glutPostRedisplay();
}

void reparam_cb(int id)
{
	float per_u = 1.0/arclen_u;
	float per_w = 1.0/arclen_w;
	float del_u = per_u * len_u;
	float del_w = per_w * len_w;
	float end_u = start_u+del_u;
	float end_w = start_w+del_w;
	if(end_u>1.0)
		end_u = 1.0;
	if(end_w>1.0)
		end_w = 1.0;
	if(surface_type==0)
	{
		generate_bezier_surface(start_u,end_u,start_w,end_w);
		cur_banner_id = 34;
	}
	glutPostRedisplay();
}

void set_ui_elems(GLUI* panel_id)
{
	GLUI_Panel *ct_panel = new GLUI_Panel( panel_id, "Curve Type" );
	GLUI_RadioGroup *radio_button_group = panel_id->add_radiogroup_to_panel(ct_panel,&curve_type,-1,radio_button_cb_func);
	panel_id->add_radiobutton_to_group(radio_button_group,"Bezier");
	panel_id->add_radiobutton_to_group(radio_button_group,"Cubic B-Spline");
	panel_id->add_radiobutton_to_group(radio_button_group,"Sub-division");
	panel_id->add_radiobutton_to_group(radio_button_group,"Rational Bezier");
	panel_id->add_radiobutton_to_group(radio_button_group,"Arbit B-Spline");
	panel_id->add_radiobutton_to_group(radio_button_group,"NURBS");

	GLUI_Panel *bc_panel = new GLUI_Panel(panel_id,"Bezier OPTS");
	GLUI_Spinner *num_samples_spinner = panel_id->add_spinner_to_panel(bc_panel,"Samples",GLUI_SPINNER_INT,&num_subdivisions,-1,samples_cb_func);//,NULL);
	num_samples_spinner->set_int_limits(10,50,GLUI_LIMIT_WRAP);
	panel_id->add_button_to_panel(bc_panel,"C1-Split",0,c1_split_cb);
	panel_id->add_button_to_panel(bc_panel,"Closed",0,bezier_close_cb);



	GLUI_Panel *sdc_panel = new GLUI_Panel(panel_id,"Sub-div OPTS");
	GLUI_Spinner *numpnt_spinner = panel_id->add_spinner_to_panel(sdc_panel,"Subdiv. Iters",GLUI_SPINNER_INT, &subdiv_iters,-1,spinner_cb_func);
 	numpnt_spinner->set_int_limits( 1, 8, GLUI_LIMIT_WRAP );

 	GLUI_Panel *rbc_panel = new GLUI_Panel(panel_id,"Rat-Bez OPTS");
 	panel_id->add_button_to_panel(rbc_panel,"Load Weights",0,load_rbsweight_cb_func);

 	GLUI_Panel *arbbsp_panel = new GLUI_Panel(panel_id,"Arb. BS OPTS");
 	GLUI_Spinner *arbbsp_spinner = panel_id->add_spinner_to_panel(arbbsp_panel,"K",GLUI_SPINNER_INT, &arbbsp_k);//,-1,spinner_cb_func);
 	arbbsp_spinner->set_int_limits( 1, 8, GLUI_LIMIT_WRAP );
 	panel_id->add_button_to_panel(arbbsp_panel,"Load Knot Vector",0,load_knotvectors_cb_func);
 	panel_id->add_button_to_panel(arbbsp_panel,"Gen Curve",0,gen_arbbsp_cb);


 	GLUI_Panel *nurb_panel = new GLUI_Panel(panel_id,"NURBS OPTS");
 	GLUI_Spinner *nurbdeg_spinner = panel_id->add_spinner_to_panel(nurb_panel,"K",GLUI_SPINNER_INT, &nurb_k);//,-1,spinner_cb_func);
 	nurbdeg_spinner->set_int_limits( 1, 8, GLUI_LIMIT_WRAP );
 	panel_id->add_button_to_panel(nurb_panel,"Load Weights",0,load_nurbweight_cb_func);
 	panel_id->add_button_to_panel(nurb_panel,"Gen Curve",0,gen_nurbs_cb);

	GLUI_Panel *gen_opt_panel = new GLUI_Panel(panel_id,"Gen OPTS");
	GLUI_RadioGroup *select_2daxis_radio = panel_id->add_radiogroup_to_panel(gen_opt_panel,&twodim_axis,-1,twodim_axis_radio_cb_func);
	panel_id->add_radiobutton_to_group(select_2daxis_radio,"X-Y Axis");
	panel_id->add_radiobutton_to_group(select_2daxis_radio,"Y-Z Axis");
	panel_id->add_checkbox_to_panel(gen_opt_panel,"Local Control Enable",&local_control_enabled);
	panel_id->add_checkbox_to_panel(gen_opt_panel,"Gridlines",&enable_grid_lines,-1,simple_redisplay_cb_func);
	panel_id->add_checkbox_to_panel(gen_opt_panel,"Draw Control Polygon",&enable_control_polygon,-1,simple_redisplay_cb_func);
	panel_id->add_button_to_panel(gen_opt_panel,"Reset",0,reset_cb_func);
	panel_id->add_button_to_panel(gen_opt_panel,"Quit", 0,(GLUI_Update_CB)exit);

	panel_id->add_column(true);
	GLUI_Panel *sor_panel = new GLUI_Panel(panel_id,"SOR OPTS");
	panel_id->add_checkbox_to_panel(sor_panel,"X-Axis",&use_x_axis);
	GLUI_Spinner *sorslice_spinner = panel_id->add_spinner_to_panel(sor_panel,"Slices",GLUI_SPINNER_INT,&num_slices);
	sorslice_spinner->set_int_limits(10,20,GLUI_LIMIT_WRAP);
	panel_id->add_button_to_panel(sor_panel,"Generate Solid",0,generate_sor_solid);

	GLUI_Panel* soe_panel = new GLUI_Panel(panel_id,"SOE OPTS");
	GLUI_Spinner *soeslice_spinner =  panel_id->add_spinner_to_panel(soe_panel,"Slices",GLUI_SPINNER_INT,&num_slices_soe);
	soeslice_spinner->set_int_limits(10,20,GLUI_LIMIT_WRAP);
	panel_id->add_edittext_to_panel(soe_panel,"Depth",GLUI_EDITTEXT_FLOAT,&depth);
	panel_id->add_edittext_to_panel(soe_panel,"x",GLUI_EDITTEXT_FLOAT,&x_dir);
	panel_id->add_edittext_to_panel(soe_panel,"y",GLUI_EDITTEXT_FLOAT,&y_dir);
	panel_id->add_edittext_to_panel(soe_panel,"z",GLUI_EDITTEXT_FLOAT,&z_dir);
	panel_id->add_button_to_panel(soe_panel,"Generate Solid",0,generate_soe);
	//panel_id->add

	GLUI_Panel* sowe_panel = new GLUI_Panel(panel_id,"SWEEP OPTS");
	panel_id->add_button_to_panel(sowe_panel,"Generate Solid",0,generate_sweep_solid);

	GLUI_Panel* fileio_panel = new GLUI_Panel(panel_id,"FILE IO");
	panel_id->add_checkbox_to_panel(fileio_panel,"Compatibility",&is_compatible);
	panel_id->add_button_to_panel(fileio_panel,"Gen. File",0,generate_file);

	GLUI_Panel* loft_panel = new GLUI_Panel(panel_id,"Lofting");
	panel_id->add_button_to_panel(loft_panel,"Draw Next",0,switch_drawing_context);
	panel_id->add_button_to_panel(loft_panel,"Generate Loft",0,generate_loft_surface);

	panel_id->add_column(true);

	GLUI_Panel* surface_panel = new GLUI_Panel(panel_id,"Surfaces");
	GLUI_RadioGroup *select_surface_type = panel_id->add_radiogroup_to_panel(surface_panel,&surface_type);//,-1,surface_selection_cb);
	panel_id->add_radiobutton_to_group(select_surface_type,"Bezier Surface");
	panel_id->add_radiobutton_to_group(select_surface_type,"B-Spline Surface");
	panel_id->add_radiobutton_to_group(select_surface_type,"T-Nurbs Surface");
	panel_id->add_button_to_panel(surface_panel,"Generate",0,generate_surface_cb);

	GLUI_Panel* bez_surface_opts_panel = new GLUI_Panel(panel_id,"Bez. S OPTS");
	GLUI_Spinner *num_samples_u_spinner = panel_id->add_spinner_to_panel(bez_surface_opts_panel,"Samples-U",GLUI_SPINNER_INT,&num_samples_in_u,-1,bezier_surface_sample_change_cb);//NULL);
	num_samples_u_spinner->set_int_limits(11,50,GLUI_LIMIT_WRAP);
	GLUI_Spinner *num_samples_w_spinner = panel_id->add_spinner_to_panel(bez_surface_opts_panel,"Samples-W",GLUI_SPINNER_INT,&num_samples_in_v,-1,bezier_surface_sample_change_cb);//NULL);
	num_samples_w_spinner->set_int_limits(11,50,GLUI_LIMIT_WRAP);
	panel_id->add_checkbox_to_panel(bez_surface_opts_panel,"G-1 Cont.",&enforce_g1_cont);

	GLUI_Panel *bsp_surface_opts_panel = new GLUI_Panel(panel_id,"BSP OPTS");
	panel_id->add_checkbox_to_panel(bsp_surface_opts_panel,"Closed in M",&draw_closed_bsp_surface_m);
	panel_id->add_checkbox_to_panel(bsp_surface_opts_panel,"Closed in N",&draw_closed_bsp_surface_n);

	GLUI_Panel *tnurbs_surface_opts_panel = new GLUI_Panel(panel_id,"T-NURBS OPTS");
	panel_id->add_edittext_to_panel(tnurbs_surface_opts_panel,"K ",GLUI_EDITTEXT_INT,&k_s);
	panel_id->add_edittext_to_panel(tnurbs_surface_opts_panel,"L ",GLUI_EDITTEXT_INT,&l_s);
	panel_id->add_checkbox_to_panel(tnurbs_surface_opts_panel,"Random Weight",&gen_random_weight_tnurbs);

	GLUI_Panel * subd_opts_panel = new GLUI_Panel(panel_id,"Subdivisions");
	GLUI_RadioGroup *select_subdiv_type = panel_id->add_radiogroup_to_panel(subd_opts_panel,&subdiv_type);
	panel_id->add_radiobutton_to_group(select_subdiv_type,"Doo-Sabine");
	panel_id->add_radiobutton_to_group(select_subdiv_type,"Catmull-Clark");
	panel_id->add_radiobutton_to_group(select_subdiv_type,"Loop");
	panel_id->add_button_to_panel(subd_opts_panel,"Subdivide",0,subdivision_processing_cb);


	GLUI_Panel* surface_gen_opts_panel = new GLUI_Panel(panel_id,"Gen OPTS");
	panel_id->add_checkbox_to_panel(surface_gen_opts_panel,"Control Polyhedron",&enable_control_polyhedron,-1,simple_redisplay_cb_func);
	panel_id->add_button_to_panel(surface_gen_opts_panel,"Read OFF",0,read_off_cb);
	panel_id->add_button_to_panel(surface_gen_opts_panel,"Export OBJ",0,export_shape_obj_cb);
	panel_id->add_edittext_to_panel(surface_gen_opts_panel,"Arc Len u",GLUI_EDITTEXT_FLOAT,&arclen_u);
	panel_id->add_edittext_to_panel(surface_gen_opts_panel,"Arc Len w",GLUI_EDITTEXT_FLOAT,&arclen_w);
	panel_id->add_edittext_to_panel(surface_gen_opts_panel,"Start u",GLUI_EDITTEXT_FLOAT,&start_u);
	panel_id->add_edittext_to_panel(surface_gen_opts_panel,"Start w",GLUI_EDITTEXT_FLOAT,&start_w);
	panel_id->add_edittext_to_panel(surface_gen_opts_panel,"Len u",GLUI_EDITTEXT_FLOAT,&len_u);
	panel_id->add_edittext_to_panel(surface_gen_opts_panel,"Len w",GLUI_EDITTEXT_FLOAT,&len_w);
	panel_id->add_button_to_panel(surface_gen_opts_panel,"Reparam.",0,reparam_cb);

	panel_id->add_column(true);

	GLUI_Panel *io_opts_panel = new GLUI_Panel(panel_id,"Input Mode");
	GLUI_RadioGroup *select_io_mode = panel_id->add_radiogroup_to_panel(io_opts_panel,&io_type);
	panel_id->add_radiobutton_to_group(select_io_mode,"Free-Form");
	panel_id->add_radiobutton_to_group(select_io_mode,"Reconstruction");

	GLUI_Panel *recons_opts_panel = new GLUI_Panel(panel_id,"Recons. OPT");
	GLUI_RadioGroup *recons_sel_mode = panel_id->add_radiogroup_to_panel(recons_opts_panel,&recons_type);
	panel_id->add_radiobutton_to_group(recons_sel_mode,"Crust");
	panel_id->add_radiobutton_to_group(recons_sel_mode,"NN-Crust");
	panel_id->add_checkbox_to_panel(recons_opts_panel,"Always OPEN",&always_open_recons);
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_RGBA|GLUT_DEPTH|GLUT_DOUBLE);
	glutInitWindowSize(1610,800);
	int window_id = glutCreateWindow("CSE5543-Lab2");

	glewExperimental = GL_TRUE;
	glewInit();

	glutDisplayFunc(display);
	glutMotionFunc(motion_func);
	glutPassiveMotionFunc(passive_motion_func);
	glutKeyboardFunc(keybd_func);
	glutMouseFunc(mouse_func);

	GLUI_Master.set_glutReshapeFunc(NULL);
	GLUI_Master.set_glutIdleFunc(NULL);

	glui_window = GLUI_Master.create_glui( "Control Panel" );
 	set_ui_elems(glui_window);
 	glui_window->set_main_gfx_window( window_id );

 	init();

	glutMainLoop();

	return 0;
}
