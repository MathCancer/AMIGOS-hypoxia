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

#include "../BioFVM.h" 


using namespace BioFVM;

int omp_num_threads = 8; // set number of threads for parallel computing
// set this to # of CPU cores x 2 (for hyperthreading)

double max_substrate_1=38.0;
void supply_function( Microenvironment* microenvironment, int voxel_index , std::vector<double>* write_here )
{
	double min_x=microenvironment->mesh.bounding_box[0];
	double max_x=microenvironment->mesh.bounding_box[3];
	double min_y=microenvironment->mesh.bounding_box[1];
	double max_y=microenvironment->mesh.bounding_box[4];
	double min_z=microenvironment->mesh.bounding_box[2];
	double max_z=microenvironment->mesh.bounding_box[5];
	double strip_width=40;
	(*write_here)[0] = 0; 

	if( abs(max_x-microenvironment->voxels(voxel_index).center[0]) < strip_width || abs(microenvironment->voxels(voxel_index).center[0]- min_x)< strip_width  
			|| abs(max_y-microenvironment->voxels(voxel_index).center[1]) < strip_width || abs(microenvironment->voxels(voxel_index).center[1]- min_y)< strip_width  
				|| abs(max_z-microenvironment->voxels(voxel_index).center[2]) < strip_width || abs(microenvironment->voxels(voxel_index).center[2]- min_z)< strip_width )
				{
					(*write_here)[0] = max_substrate_1;
				}	
				
	return ;
}

void supply_target_function( Microenvironment* m, int voxel_index , std::vector<double>* write_here )
{
	(*write_here)[0] = max_substrate_1; 	
	return ;
}


int writeCellReport(std::vector<Basic_Agent*> all_cells, double timepoint, double radius)
{
	std::string filename; 
	filename.resize( 1024 ); 
	sprintf( (char*) filename.c_str() , "cells_%i_conv_test4_1.txt" , (int)timepoint ); 
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	povFile<<"\tID\tx\ty\tz\tradius\tphenotype\n";
	int phenotype_code;
	for(int i=0;i<all_cells.size();i++)
	{
		phenotype_code=-1;
		povFile<<i<<"\t"<<all_cells[i]->ID<<"\t"<<all_cells[i]->position[0]<<"\t" << all_cells[i]->position[1] <<"\t"<< all_cells[i]->position[2]<<"\t";
		povFile<<radius<<"\t"<<phenotype_code<<std::endl;
	}
	povFile.close();
	return 0;
}


void process_output(double t, double dt, double mesh_resolution, Microenvironment microenvironment)
{
	std::cout << "current simulated time: " << t   << " minutes " << std::endl; 
	std::cout << "interval wall time: ";
	BioFVM::TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::stopwatch_value() ); 
	std::cout << std::endl; 
	std::cout << "total wall time: "; 
	BioFVM::RUNTIME_TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
	std::cout << std::endl;
	
	std::cout << "time: "<<t<<std::endl;
	BioFVM::TIC();
	
	std::string filename; 
	filename.resize( 1024 , '\0' ); 
	sprintf( (char*) filename.c_str() , "output_%f_%f_%f.mat" , t, dt, mesh_resolution ); 
	filename.resize( strlen( filename.c_str() ) ); 
	std::cout << "\tWriting to file " << filename << " ... " << std::endl; 
	microenvironment.write_to_matlab( filename ); 
}

int main( int argc, char* argv[] )
{
	double t = 0.0; 
	double t_output_interval = 0.5; 
	double t_next_output_time = 0;
	int next_output_index = 0;
	double tumor_focus_radius = 150;
	
	
	double dt = strtod(argv[1],NULL);
	double t_max = strtod(argv[2],NULL);
	double mesh_resolution=strtod(argv[3],NULL);
	double cell_uptake_rate=10.0;
	
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.set_density(0, "substrate1" , "dimensionless" );
	
	std::vector<double> center(3);
	center[0] = 0; center[1] = 0; center[2] = 0; 
	
	double minX=-500, minY=-500, minZ=-500, maxX=500, maxY=500, maxZ=500;//, mesh_resolution=10; 
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	microenvironment.diffusion_coefficients[0] = 100000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.1;
	microenvironment.display_information( std::cout );

	microenvironment.bulk_supply_rate_function =  supply_function;
	microenvironment.bulk_supply_target_densities_function = supply_target_function;
	
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.density_vector(i)[0]= max_substrate_1; 
	}
	
	
	for(int i=0;i<microenvironment.number_of_voxels();i++)
	{
		if(sqrt(norm_squared(microenvironment.voxels(i).center))<tumor_focus_radius)
		{
			Basic_Agent * temp_point_sink = create_basic_agent();
			temp_point_sink->register_microenvironment(&microenvironment);
			temp_point_sink->assign_position(microenvironment.voxels(i).center);
			temp_point_sink->set_total_volume( pow(mesh_resolution,3.0) );
			(*temp_point_sink->uptake_rates)[0]=cell_uptake_rate;
			temp_point_sink->set_internal_uptake_constants(dt); 
		}
	}
	
	writeCellReport(all_basic_agents,0.0, 10.0);
	// display information 
	microenvironment.display_information( std::cout );
	int center_voxel_index = microenvironment.nearest_voxel_index(center);
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	int output_index=0;	
	
	std::cout<<"number of agents: "<< all_basic_agents.size() <<std::endl;
	
	while( t < t_max )
	{		
		if(  fabs( t - t_next_output_time ) < dt/10.0 )
		{
			process_output(t, dt, mesh_resolution, microenvironment);
			t_next_output_time += t_output_interval; 
			next_output_index++;  
		}
		// simulate microenvironment 
		microenvironment.simulate_bulk_sources_and_sinks( dt );
		microenvironment.simulate_cell_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		
		t += dt; 
		output_index++; 
	}
	process_output(t_max, dt, mesh_resolution, microenvironment);
	std::cout << "done!" << std::endl; 
	return 0; 
}