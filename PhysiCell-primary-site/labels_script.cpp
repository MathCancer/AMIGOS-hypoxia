/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.3.1) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.3.1) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

// set number of threads for OpenMP (parallel computing)
int omp_num_threads = 4; // set this to # of CPU cores x 2 (for hyperthreading)

#include "./BioFVM/BioFVM_matlab.h" 

void minutes_to_label( char* str, double t ); 

int main( int argc, char* argv[] )
{
	// OpenMP setup
	omp_set_num_threads(omp_num_threads);
	
	char base_command [1024];
	strcpy( base_command , "magick mogrify -font Arial -fill black -pointsize 75 -gravity NorthWest -annotate +20+20" ); 
		
	#pragma omp parallel for 
	for( int i= atoi( argv[2] ) ; i <= atoi( argv[3] ); i++ )
	{
		char str [1024]; 
		char png_filename [1024];
		char mat_filename [1024]; 
		
		char label [2048]; 

		sprintf( png_filename , "pov%08u.png" , i ); 
		sprintf( mat_filename , "%s/output%08i_cells_physicell.mat" , argv[1] , i );
		
		std::vector< std::vector<double> > MAT = BioFVM::read_matlab( mat_filename );
		
		int number_of_cells = MAT[0].size(); 
		double time_in_minutes = (double) i * 60.0; 
		
		minutes_to_label( str, time_in_minutes ); 
		/*
		std::cout << str << std::endl; 
		std::cout << number_of_cells << std::endl; 
		std::cout << mat_filename << std::endl; 
		std::cout << png_filename << std::endl << std::endl; 
		*/
		
		sprintf( label, "%s\\n%i agents" , str , number_of_cells ) ; 
		
//		std::cout << label << std::endl << std::endl; 
		
		char my_command [2048]; 
		sprintf( my_command, "%s \"%s\" %s" , base_command , label , png_filename ); 
		
		std::cout << my_command << std::endl << std::endl; 
		
		system( my_command ); 	
	}
	

	
	return 0; 
}

void minutes_to_label( char* str, double t )
{
	int days = floor( t / (24 * 60 )); 
	t -= 24*60*days; 

	int hours = floor( t/ 60); 
	t -= hours*60;  

	double minutes = t; 

	sprintf( str, "Current time: %i days, %i hours, and %3.2f minutes" , days, hours, minutes );
	return;
}
