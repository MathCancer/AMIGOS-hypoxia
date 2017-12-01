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
double pi= 3.1415926535897932384626433832795;
int main( int argc, char* argv[] )
{
	double t = 0.0; 
	double dt = strtod(argv[1],NULL);
	double t_max = strtod(argv[2],NULL);
	double mesh_resolution=strtod(argv[3],NULL);
	
	std::cout<<"dt: "<< dt<<", tMax: "<< t_max<<", dx:"<< mesh_resolution<<std::endl;
	// openmp setup
	omp_set_num_threads(omp_num_threads);
	
	// create a microenvironment; 
	Microenvironment microenvironment;
	microenvironment.set_density(0, "substrate1" , "dimensionless" );
	
	std::vector<double> center(3);
	center[0] = 0; center[1] = 0; center[2] = 0; 
	
	double minX=-500, minY=-mesh_resolution/2, minZ=-mesh_resolution/2, maxX=500, maxY=mesh_resolution/2, maxZ=mesh_resolution/2; 
	microenvironment.resize_space_uniform( minX,maxX,minY,maxY,minZ,maxZ, mesh_resolution );
	// register the diffusion solver 	
	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D; 
	
	microenvironment.diffusion_coefficients[0] = 100000; // microns^2 / min 
	microenvironment.decay_rates[0] = 0.0;
	
	
				
	// display information 
	microenvironment.display_information( std::cout );
	
	int center_voxel_index = microenvironment.nearest_voxel_index(center);
	
	t=0;
	
	double L0=500.0, D=microenvironment.diffusion_coefficients[0];
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		double new_value;
		double X = sqrt(norm_squared(microenvironment.voxels(i).center ));
		new_value=1.0+ cos((pi/L0) *X)* exp(-(D* pi*pi/ (L0*L0)) *t);
		microenvironment.density_vector(i)[0]=new_value;
	}
	
	while( t < t_max )
	{
		microenvironment.simulate_diffusion_decay( dt );
		t += dt; 
	}
	microenvironment.write_to_matlab( "final.mat" );
	std::cout << "done!" << std::endl; 
	return 0; 
}