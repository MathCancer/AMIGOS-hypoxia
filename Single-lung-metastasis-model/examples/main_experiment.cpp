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

//setup Mersenne Twister random generator 
std::random_device rd;
std::mt19937 gen(rd());

const int TUMOR_TYPE= 1;
const int VESSEL_TYPE= 2;

const double hypoxic_thresh1=8;
const double hypoxic_thresh2=15;

const double max_vegf_secretion_rate=1;

double UniformRandom()
{
	return std::generate_canonical<double, 10>(gen);
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
	double dt = 0.01; 
	double t_output_interval = 60.0; // 1.0; 
	double t_max = 30.0 * 24.0 * 60.0; // 500.0;  
	double t_next_output_time = 0; 
	int next_output_index = 0; 
	
	double dx; 
	double dy;  
	double dz; 
	
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// set random seed 
	srand(4); 
	
	// figure out the bounding box 
	std::vector<double> bounding_box( 6, 0.0 );
	// min_x_index=0, min_y_index=1, min_z_index=2, max_x_index=3,	max_y_index=4, max_z_index=5
	bounding_box[0] = 0; bounding_box[3] = 5000; 
	bounding_box[1] = 0; bounding_box[4] = 5000; 
	bounding_box[2] = 0; bounding_box[5] = 5000; 
	dx=20; dy=dx; dz=dx;

	
	// create a Microenvironment; 
	BioFVM::Microenvironment microenvironment;
	microenvironment.name="substrate scale";
	microenvironment.set_density(0, "oxygen" , "mmHg" );
	microenvironment.add_density( "vegF" , "dimensionless" );
	
	std::cout << bounding_box << std::endl; 
	
	microenvironment.resize_space( bounding_box[0] , bounding_box[3] , bounding_box[1], bounding_box[4] , bounding_box[2] , bounding_box[5] ,dx,dy,dz );
	microenvironment.spatial_units = "microns";
	microenvironment.time_units = "minutes";
	microenvironment.mesh.units = "microns";	

	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		microenvironment.density_vector(n)[0] = 38.0; 	
		microenvironment.density_vector(n)[1] = 0.0; 
	}
	
	// register the diffusion solver 
	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
		
	// register bulk properties 
	microenvironment.diffusion_coefficients[0] = 100000; // microns^2 / min 
	microenvironment.diffusion_coefficients[1] = 1000 ; 

	microenvironment.decay_rates[0] = 0.1;
	microenvironment.decay_rates[1] = 0.01;
	
	std::cout<< "Decay rates: "<< microenvironment.decay_rates[0]<<", "<<microenvironment.decay_rates[1]<<std::endl;
	
	// display information 	
	microenvironment.display_information( std::cout );
	
	std::vector< std::vector<double> > tumor_cells_position = read_matlab( "tumor_cells.mat" );
	double cell_radius=10;
	for(int i=0; i< tumor_cells_position.size();i++)
	{
		std::vector<double> tempPoint(3,0.0);
		for( int j=0; j < 3 ; j++ )
		{ tempPoint[j] = tumor_cells_position[i][j]; }		

		Basic_Agent* cancer_cell = create_basic_agent();
		cancer_cell->register_microenvironment(&microenvironment);
		cancer_cell->assign_position(tempPoint);
		cancer_cell->set_total_volume( (4.0/3.0)*pi*pow(cell_radius,3.0) );
		(*cancer_cell->uptake_rates)[0]=10;
		cancer_cell->set_internal_uptake_constants(dt); 
		cancer_cell->type= TUMOR_TYPE;
	}

		
	std::vector< std::vector<double> > vessels_cells_position = read_matlab( "blood_cells.mat" );
	
	for(int i=0; i< vessels_cells_position.size();i++)
	{
		std::vector<double> tempPoint(3,0.0);
		for( int j=0; j < 3 ; j++ )
		{ tempPoint[j] = vessels_cells_position[i][j]; }		
		Basic_Agent* vessel_cell = create_basic_agent();
		vessel_cell->register_microenvironment(&microenvironment);
		vessel_cell->assign_position(tempPoint);
		vessel_cell->set_total_volume( (4.0/3.0)*pi*pow(cell_radius,3.0) );
		(*vessel_cell->secretion_rates)[0]=10;
		(*vessel_cell->saturation_densities)[0]=70;
		vessel_cell->set_internal_uptake_constants(dt); 
		vessel_cell->type= VESSEL_TYPE;
	}
	
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	int output_index=0;
	int num_new_cells= all_basic_agents.size();
	int num_deaths=0;
	try {
	while( t < t_max )
	{
		#pragma omp parallel for
		for(int i=0;i<all_basic_agents.size();i++){
			if(all_basic_agents[i]->type==TUMOR_TYPE)
			{
				(*all_basic_agents[i]->secretion_rates)[1] = 0;
				double o2_conc=(all_basic_agents[i]->nearest_density_vector())[0] ;
				if(o2_conc< hypoxic_thresh2)
				{
					if(o2_conc<=hypoxic_thresh1)
						(*all_basic_agents[i]->secretion_rates)[1] = max_vegf_secretion_rate;
					else
						(*all_basic_agents[i]->secretion_rates)[1] = ((hypoxic_thresh2- o2_conc)/(hypoxic_thresh2-hypoxic_thresh1))*max_vegf_secretion_rate;
					(*all_basic_agents[i]->saturation_densities)[1]=1;
				}
				all_basic_agents[i]->set_internal_uptake_constants(dt); 
			}
		}
	
	
		if(  fabs( t - t_next_output_time ) < 0.001 )
		{
			process_output(t, dt, dx, microenvironment);
			t_next_output_time += t_output_interval; 
			next_output_index++; 
		}
		// simulate microenvironment 
		microenvironment.simulate_cell_sources_and_sinks( dt );
		microenvironment.simulate_diffusion_decay( dt );
		t += dt; 
		output_index++; 
	}
	process_output(t_max, dt, dx, microenvironment);
	BioFVM::RUNTIME_TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );
	std::cout << "done!" << std::endl; 
	}
	catch( const std::exception& e ) { 
		std::cout << e.what(); 
	}
	return 0; 
}
