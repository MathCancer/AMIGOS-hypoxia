/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.5) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#ifndef __BioFVM_mesh_h__
#define __BioFVM_mesh_h__

#include <iostream>
#include <vector> 

#include "BioFVM_matlab.h"

namespace BioFVM{

 /*! \brief Voxels are the basic spatial container for densities, which are networked into meshes. 
  * 
  * Voxels are the basic spatial container for a finite volume method. Voxels are connected to 
  * other voxels into a General_Mesh, here most likely a Cartesian_Mesh. Voxel boundaries are ?.
  * 
  * A Microenvironment Domain will include a network of Voxels (and ?),
  * with a vector<double> of densities for each Voxel, along with rate constants, etc. 
  * The Domain may also include a vector<double> of flux coefficients for each ?.
 */
    
// How long have these notes been in your code PM?
 
class Voxel
{

 private:
	friend std::ostream& operator<<(std::ostream& os, const Voxel& v); 
	/*!< outputs the Voxel to an open ostream 
	 * \param os -- the stream 
	 * \param mv -- the voxel you use this friendly friend operator on
	 * Example: Voxel v; 
	 *          cout << v << endl; 
	*/ 

 public:
	Voxel(); 
	int mesh_index; /*!< voxel's index in a General_Mesh */ 

	double volume; /*!< voxel's volume (cubic spatial units) */ 
	std::vector<double> center; /*!< center of volume */  // x, y, z
	bool is_Dirichlet;
	void stream_output_with_units( std::ostream& os , std::string units ) const;
    
};

class Voxel_Link
 {
    public:
    int ID;
    double surface_area;
    Voxel* pVoxel1;
    Voxel* pVoxel2;
    double center_to_center_distance;
     std::vector<double> normal_vector_i_to_j;
    Voxel_Link();
//     ~Voxel_Link();
    Voxel_Link(Voxel* pVoxel1, Voxel* pVoxel2, double surface_area, double center_to_center_distance,  std::vector<double>  normal_vector_i_to_j, int ID);
 };
    
// NOT USING FACES
    
//class Voxel_Face
//{
// private:
//	friend std::ostream& operator<<(std::ostream& os , const Voxel_Face& vf ); 
//	
// public:
//	Voxel_Face(); 
//	int mesh_index; 
//	
//	double surface_area; 
//	std::vector<double> center; 
//	std::vector<double> outward_normal; 
//	std::vector<double> inward_normal; 
//	
//	void stream_output_with_units( std::ostream& os , std::string units ) const;
//};

class General_Mesh
{
 private: 
	friend std::ostream& operator<<(std::ostream& os, const General_Mesh& mesh);  
	
	// this stores the indexing of the voxel faces (connect voxel i to voxel j, face stored at k)
	// only for use in a future release
	// std::unordered_map< int,std::unordered_map<int,int> > voxel_face_index_mapping; 
	
 public:
	General_Mesh();  
	
    void activate_voxel_links( void );  // Done
    void initialize_connected_voxel_indices (void);  // Should this make the voxel links too?
    int add_voxel (double x_center, double y_center, double z_center, double volume); // Done.  Assigns mesh index based on the size of voxels.  [0] of the voxels vector will be the first voxel created.  That voxel will be index 0.
    int max_voxel_to_voxel_connections;

    int heaviside_fn ( double in );
    
	// [xmin ymin zmin xmax ymax zmax ]
	std::vector<double> bounding_box; 
	
	std::vector<Voxel> voxels;
    std::vector<Voxel_Link> voxel_links; // Vector of voxel links - should this be a double or a set of points or a vector?
    //std::vector<double> center_to_center_voxel_distances; // Holds center to center distances between all linked voxels.  // NOT IN USE - CURRENTLY IN TISSUE ONLY, COULD BE USEUFL HERE ALSO, BUT MOVING ON! CONSIDER PUTTING THE MEASURE STRAIGHT INTO THE VOXEL LINKS
	// each voxel[k] has a list of connected voxels -- helpful for some numerical methods
	std::vector< std::vector<int> > connected_voxel_indices;
	
	int nearest_voxel_index( std::vector<double>& position );   
	bool is_position_valid(double x, double y, double z);

	void connect_voxels_indices_only(int i,int j, double SA);
	
    void link_voxels(Voxel* pVoxel1, Voxel* pVoxel2);//Done
    void link_voxels(Voxel* pVoxel1, Voxel* pVoxel2, double surface_area, double center_to_center_distance, std::vector<double> normal_vector_i_to_j);
    void link_voxels(int m, int n, double surface_area, double center_to_center_distance, std::vector<double> normal_vector_i_to_j); // Done
    
	/*! This removes all connections between voxels[i] and voxels[j], and deletes the associated 
	    Voxel_Face(s). */
    
    // Do I need to have a disconnector for the links also?
    
    
//	void disconnect_voxels(int i, int j);
//	void clear_voxel_face_index_mapping( void );
	
	bool Cartesian_mesh; 
	bool uniform_mesh; 
	bool regular_mesh;
	bool use_voxel_links;
	
	std::string units; 
	
	void display_information( std::ostream& os); 
	
	void write_to_matlab( std::string filename );
    void write_mesh_to_matlab( std::string filename );  // Done
    void write_links_to_matlab( std::string filename );  // Done

	void read_from_matlab( std::string filename ); 
};

class Cartesian_Mesh : public General_Mesh
{
 private:
 
 public:
	std::vector<double> x_coordinates; 
	std::vector<double> y_coordinates;
	std::vector<double> z_coordinates;
	std::vector< std::vector<int> > moore_connected_voxel_indices; // Keeps the list of voxels in the Moore nighborhood
	void create_moore_neighborhood(void);
	int voxel_index( int i, int j, int k ); 
	std::vector<int> cartesian_indices( int n ); 
	
	double dx;
	double dy;
	double dz; 
	
	double dV; 	
	double dS;

	double dS_xy;
	double dS_yz; 
	double dS_xz;
	
	Cartesian_Mesh(); // done 
	
	Cartesian_Mesh( int , int , int );  
    void link_voxels(int m, int n, double surface_area);
    bool use_voxel_links;
    // delete these notes
/*
    in the constructors, and resize operations
    
    if( use_voxel_links == true )
    {
        create_voxel_links();
    }
 */   
    // end of that comment
    
    void create_voxel_links( void );  // done
    
	void resize( int,int,int ); 
	void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , int x_nodes, int y_nodes, int z_nodes ); 
	void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz ); 
	void resize_uniform( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx ); 
	
	int nearest_voxel_index( std::vector<double>& position );   
//	int nearest_voxel_face_index( std::vector<double>& position );
	std::vector<int> nearest_cartesian_indices( std::vector<double>& position ); 
	Voxel& nearest_voxel( std::vector<double>& position ); 
	
	void display_information( std::ostream& os ); 
	
	void read_from_matlab( std::string filename ); 
};

class Voronoi_Mesh : public General_Mesh
{
 private:
 
 public:
	void display_information( std::ostream& os); 
};

    
};

#endif
