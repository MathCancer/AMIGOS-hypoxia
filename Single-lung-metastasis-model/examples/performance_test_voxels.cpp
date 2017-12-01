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

int num_substrates;
int center_voxel_index=0;

void supply_function( Microenvironment* microenvironment, int voxel_index , std::vector<double>* write_here )
{
	if( voxel_index==center_voxel_index)
	{	
		(*write_here)[0] = 100;
	}	
				
	return ;
}

void supply_target_function( Microenvironment* m, int voxel_index , std::vector<double>* write_here )
{
	if( voxel_index==center_voxel_index)
	{	
		(*write_here)[0] = 100;
	}
	return ;
}

int write_report(int num_voxels, double time)
{
	std::ofstream report ("scaling_test_results_voxels.txt", std::ofstream::app);
	report<<num_voxels<<"\t"<<time<<std::endl;
	report.close();
	return 0;
}

void process_output(double t, double dt, double mesh_resolution)
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
}

int main( int argc, char* argv[] )
{
	double t = 0.0; 
	double t_output_interval = 1.0; 
	double t_next_output_time = 0;
	int next_output_index = 0;
	
	double dt = 0.01;
	double t_max = 2.0; 
	double mesh_resolution=10.0;
	double domain_half_side= strtod(argv[1], NULL);
	
	num_substrates=1;
	
	std::vector<double> center(3);
	center[0] = 0; center[1] = 0; center[2] = 0; 
		 
	// openmp setup
	omp_set_num_threads(omp_num_threads);
		
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.set_density(0, "substrate0" , "dimensionless" );
	microenvironment.diffusion_coefficients[0] = 100000; 
	microenvironment.decay_rates[0] = 0.1;
	
	double minX=-domain_half_side, minY=-domain_half_side, minZ=-domain_half_side, maxX=domain_half_side, maxY=domain_half_side, maxZ=domain_half_side;//, mesh_resolution=10; 
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	microenvironment.bulk_supply_rate_function =  supply_function;
	microenvironment.bulk_supply_target_densities_function = supply_target_function;
	
	#pragma omp parallel for 
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.density_vector(i)[0]= 100; 
	}
	center_voxel_index= microenvironment.nearest_voxel_index(center);
	// display information 
	microenvironment.display_information( std::cout );
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	int output_index=0;	
		
	std::cout<<"number of agents: "<< all_basic_agents.size() <<std::endl;

	while( t < t_max )
	{		
		if(  fabs( t - t_next_output_time ) < dt/10.0 )
		{
			process_output(t, dt, mesh_resolution);
			t_next_output_time += t_output_interval; 
			next_output_index++; 
		}
		// simulate microenvironment 
		microenvironment.simulate_bulk_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );		
		t += dt; 
		output_index++; 
	}
	process_output(t_max, dt, mesh_resolution);
	BioFVM::RUNTIME_TOC();
	write_report(microenvironment.number_of_voxels(), BioFVM::runtime_stopwatch_value());
	std::cout << "done!" << std::endl; 
	return 0; 
}