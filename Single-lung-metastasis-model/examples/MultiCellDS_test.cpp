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

#include "BioFVM.h"

#include <omp.h>

using namespace BioFVM; 

void supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment: 
	// microenvironment->density_vector(voxel_index)[j]

	return; 
}

void supply_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment: 
	// microenvironment->density_vector(voxel_index)[j]

	return; 
}

void uptake_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment: 
	// microenvironment->density_vector(voxel_index)[j]

	return; 
}

int main( int argc, char* argv[] )
{
	omp_set_num_threads( 8 );
	
	std::cout << "Starting program ... " << std::endl;
	
	// create a microenvironment, and set units 
		
	Microenvironment M; 
	M.name = "microenvironment"; 
	M.time_units = "min"; 
	M.spatial_units = "micron"; 
	M.mesh.units = M.spatial_units;
	
	// set up and add all the densities you plan 

	M.set_density( 0 , "oxygen" , "mmHg" ); 
	
	// here's how you add a new substrate 
	M.add_density( "crayons" , "Megacrayola" ); 
	
	
	// set the properties of the diffusing substrates 
	
	M.diffusion_coefficients[0] = 1e5;   
	M.decay_rates[0] = 10; // 100 micron length scale 
	M.diffusion_coefficients[1] = 1e4;   
	M.decay_rates[1] = 1.0/9.0; // 300 micron length scale 
	
	// set the mesh size 
	
	double dx = 20; // 
	M.resize_space( -1000.0 , 1000.0 , -800.0, 800.0 , -500.0 , 1000.0 , dx, dx, dx );  
	
	// display summary information 
	
	M.display_information( std::cout ); 
	
	// set up metadata 

	BioFVM_metadata.program.program_name = "BioFVM MultiCellDS Test";
	BioFVM_metadata.program.program_version = "1.0";
	BioFVM_metadata.program.program_URL = "http://BioFVM.MathCancer.org";

	BioFVM_metadata.program.creator.surname = "Macklin";
	BioFVM_metadata.program.creator.given_names = "Paul";
	BioFVM_metadata.program.creator.email = "Paul.Macklin@usc.edu";
	BioFVM_metadata.program.creator.URL = "http://BioFVM.MathCancer.org"; 
	BioFVM_metadata.program.creator.organization = "University of Southern California";
	BioFVM_metadata.program.creator.department = "Center for Applied Molecular Medicine";
	BioFVM_metadata.program.creator.ORCID = "0000-0002-9925-0151";	
	
	
	BioFVM_metadata.program.citation.DOI = "10.1093/bioinformatics/btv730";	
	BioFVM_metadata.program.citation.PMID = "26656933";	
	BioFVM_metadata.program.citation.PMCID = "PMC1234567";	
	BioFVM_metadata.program.citation.text = "A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient parallelized diffusive transport solver for 3-D biological simulations, Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730.";	
	BioFVM_metadata.program.citation.notes = "notes here";	
	BioFVM_metadata.program.citation.URL = "http://dx.doi.org/10.1093/bioinformatics/btv730";	
	
	BioFVM_metadata.program.user.surname = "Kirk";
	BioFVM_metadata.program.user.given_names = "James T.";
	BioFVM_metadata.program.user.email = "Jimmy.Kirk@starfleet.mil";
	BioFVM_metadata.program.user.organization = "Starfleet";
	BioFVM_metadata.program.user.department = "U.S.S. Enterprise (NCC 1701)";
	BioFVM_metadata.program.user.ORCID = "0000-0000-0000-0000";
	
	BioFVM_metadata.data_citation.DOI = "10.1093/bioinformatics/btv730";	
	BioFVM_metadata.data_citation.PMID = "12345678";	
	BioFVM_metadata.data_citation.PMCID = "PMC1234567";	
	BioFVM_metadata.data_citation.text = "A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient parallelized diffusive transport solver for 3-D biological simulations, Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730.";	
	BioFVM_metadata.data_citation.notes = "notes here";	
	BioFVM_metadata.data_citation.URL = "http://dx.doi.org/10.1093/bioinformatics/btv730";	

	// set BioFVM MultiCellDS options here (defaults are retained below): 
	
	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true ); 
	
	// set initial conditions 
	
	// use this syntax to create a zero vector of length 3
	// std::vector<double> zero(3,0.0); 
	
	// use this syntax for a parallelized loop over all the 
	// voxels in your mesh: 	
	
	std::vector<double> middle(3,0.0); 
	middle[0] = 42.0; 
	middle[1] = 108.0; 
	middle[2] = -13.0; 
	
	#pragma omp parallel for 
	for( int i=0 ; i < M.number_of_voxels() ; i++ )
	{
		// use this syntax to access the coordinates (as a vector) of 
		// the ith voxel; 
		// M.mesh.voxels[i].center 
		
		// use this access the jth substrate at the ith voxel
		// M.density_vector(i)[j]

		std::vector<double> displacement = M.mesh.voxels[i].center - middle; 
		double x2 = norm_squared( displacement ); 
		
		// use this access the jth substrate at the ith voxel
		
		M.density_vector(i)[0] = x2;  
		M.density_vector(i)[1] = exp( -x2 / (2*100*100 ) );  
	}
	
	// now let's make some cells
	
	srand( time(NULL)); 
	
	double length_x = M.mesh.bounding_box[3] - M.mesh.bounding_box[0];
	double length_y = M.mesh.bounding_box[4] - M.mesh.bounding_box[1];
	double length_z = M.mesh.bounding_box[5] - M.mesh.bounding_box[2];

	for( int i=0; i < 250000 ; i++ )
	{
		Basic_Agent* pBA = create_basic_agent(); 
		pBA->register_microenvironment( &M ); 
		pBA->secretion_rates->assign( M.number_of_densities() , 3.0 ); 
		pBA->uptake_rates->assign( M.number_of_densities() , 2.0 ); 
		pBA->saturation_densities->assign( M.number_of_densities() , 1.0 ); 
		
		pBA->position[0] = M.mesh.bounding_box[0] + length_x * rand() / (double) RAND_MAX; 
		pBA->position[1] = M.mesh.bounding_box[1] + length_x * rand() / (double) RAND_MAX; 
		pBA->position[2] = M.mesh.bounding_box[2] + length_x * rand() / (double) RAND_MAX; 
	}

	// save the data
	
	TIC(); 
	// M.write_to_matlab( "sample.mat" ); // barebones -- won't include cells!
	double current_simulation_time = 10.347; 
	save_BioFVM_to_MultiCellDS_xml_pugi( "sample" , M , current_simulation_time ); // MultiCellDS digital snapshot


	// use a syntax like this if you plan to save a series of files
	
	char filename_base [1024]; 
	sprintf( &filename_base[0] , "sample_%f", current_simulation_time ); 
	save_BioFVM_to_MultiCellDS_xml_pugi( filename_base , M , current_simulation_time ); 
	TOC();
	display_stopwatch_value( std::cout , stopwatch_value() ); 
	
	
	return 0; 
}
