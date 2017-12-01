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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <random>

#include "../BioFVM.h" 

using namespace BioFVM;

int omp_num_threads = 8; // set number of threads for parallel computing
// set this to # of CPU cores x 2 (for hyperthreading)
double pi= 3.1415926535897932384626433832795;

double UniformRandom()
{	
	return ((double) rand() / (RAND_MAX));
}

int main( int argc, char* argv[] )
{	
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.name="substrate scale";

	microenvironment.set_density(0, "substrate1" , "dimensionless" );
	microenvironment.spatial_units = "microns";
	microenvironment.mesh.units = "microns";
	microenvironment.time_units = "minutes";
	
	
	double minX=0, minY=0, minZ=0, maxX=1000, maxY=1000, maxZ=1000, mesh_resolution=10;
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	microenvironment.display_information( std::cout );

	
	std::vector<double> center(3);
	center[0] = (microenvironment.mesh.bounding_box[0]+microenvironment.mesh.bounding_box[3])/2;
	center[1] = (microenvironment.mesh.bounding_box[1]+microenvironment.mesh.bounding_box[4])/2;
	center[2] = (microenvironment.mesh.bounding_box[2]+microenvironment.mesh.bounding_box[5])/2;
	double stddev_squared = -100.0 * 100.0; 	
	std::vector<double> one( microenvironment.density_vector(0).size() , 1.0 ); 
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		std::vector<double> displacement = microenvironment.voxels(i).center - center; 
		double distance_squared = norm_squared( displacement );
		double coeff = distance_squared;
		coeff /=  stddev_squared;
		microenvironment.density_vector(i)[0]= exp( coeff ); 
	}
	microenvironment.write_to_matlab( "initial_concentration.mat" ); 

	
	

	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	// register substrates properties 
	microenvironment.diffusion_coefficients[0] = 1000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.01;
				
	// display information 	
	microenvironment.display_information( std::cout );
	
	double dt = 0.01; 
	double cell_radius=5;
	for(int i=0; i< 500;i++)
	{
		std::vector<double> tempPoint(3,0.0);
		for( int j=0; j < 3 ; j++ )
		{ tempPoint[j] = UniformRandom()*1000; }		
		
		Basic_Agent * temp_point_source = create_basic_agent();
		temp_point_source->register_microenvironment(&microenvironment);
		temp_point_source->assign_position(tempPoint);
		temp_point_source->set_total_volume( (4.0/3.0)*pi*pow(cell_radius,3.0) );
		(*temp_point_source->secretion_rates)[0]=10;
		(*temp_point_source->saturation_densities)[0]=1;
		temp_point_source->set_internal_uptake_constants(dt); 
		

		for( int j=0; j < 3 ; j++ )
		{ tempPoint[j] = UniformRandom()*1000; }		
		Basic_Agent * temp_point_sink = create_basic_agent();
		temp_point_sink->register_microenvironment(&microenvironment);
		temp_point_sink->assign_position(tempPoint);
		temp_point_sink->set_total_volume( (4.0/3.0)*pi*pow(cell_radius,3.0) );
		(*temp_point_sink->uptake_rates)[0]=0.8;
		temp_point_sink->set_internal_uptake_constants(dt); 
	}
	
	
	double t = 0.0; 
	double t_max=5;

	while( t < t_max )
	{
		microenvironment.simulate_cell_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		t += dt; 
	}
	microenvironment.write_to_matlab( "final.mat" );
	std::cout<<"done!"<<std::endl;
	return 0; 
}