/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.4) [1]        #
#                                                                           #
# [1] A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient  #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2016, Paul Macklin and the BioFVM Project              #
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
//#include "PhysiCell_SVG.h"

using namespace BioFVM; // Change?

// Not currently using

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

// Not currently using the above

int main( int argc, char* argv[] )
{
	RUNTIME_TIC();
	omp_set_num_threads( 4 );
	
	std::cout << "Starting program ... " << std::endl;
    
    Phenotype breast_cancer;
    Phenotype vasculature;
    

    
    // Setting up tissue
    
    // Perhaps the voxel constructor shoudl just be overloaded to create "tissue" voxels that use the population vectors.  Why not do that?
    
   Tissue tissue;
    double cell_radius = 6.5; // um
    double cell_volume = 4.0/3.0* 3.141592654 * pow(cell_radius, 3);
    double voxel_radius = 20; // um - produces voxels containing 29 spherical voxels
    double test_distance = 2*voxel_radius;
    double domain_radius = 270; // um
    double exchange_surface_area = 418; // um^2, surface of sphere of radius 20/12 (number of neighbhors on a
    
    // Below three commands are used to make and link a hexagonal lattice carved into a circular shape
    
//    tissue.make_tumor_spheroid_2D( domain_radius , voxel_radius ); // Makes a 2D spheroid (circle) composed of multiple spheres (a bit like a bunch of HDSs) connected in a hexagonal lattice
//
//    // Activates the links structure
//
//    tissue.tissue_mesh.activate_voxel_links();
//
//    // Link voxels
//
//    tissue.link_tumor_spheroid( test_distance, 418 );  // same radius as the voxel initialixation.  Surface area is voxel SA divided by 12 - number of neighbors in a closest packing 3-D hex lattice.
    
    // USE BELOW LINE TO MAKE A SPHERICAL DOMAIN
    
    tissue.spherical_geometry_linker ( voxel_radius,  domain_radius );  // Actual 3-D concentric spheres.  MAKES DOMAIN, LINKERS Structures, and LINKS all at once.
    
    // END Spherical domain 
    
    std::cout<<tissue.tissue_mesh.voxels.size()<<std::endl;
    std::cout<<tissue.tissue_mesh.voxel_links.size()<<std::endl;
    
    // Create populations, add phenotypes, and let them grow ...
    
    tissue.create_population_vectors();
    
    breast_cancer.max_cells = 29.0;  // make this flexible but will do for the moment.
    breast_cancer.uptake_rate[0][0] = 10.2;  // O2 Update  - Live Cells
    breast_cancer.uptake_rate[0][1] = 0.0;  // Apop cells
    breast_cancer.uptake_rate[0][2] = 10.2;  // Necro cells
    breast_cancer.spatial_aggregation_density = 0.1; // Starting at 0.0.
//    breast_cancer.substrate_target[1][0] = 1.0;  // AF factor max secretion rate - not used in code for some reason - using the "base secretion rate" fro some reason.  
    
    vasculature.phenotype_ID = 1;
    vasculature.name = "Vasculature";
    vasculature.birth_rate = 1.0/18.0/60;
    vasculature.apoptotic_death_rate = 0.0;
    vasculature.apoptotic_clearance_rate =0.0;
    vasculature.necrotic_death_rate = 0.0;
    vasculature.necrotic_clearance_rate =0.0;
    vasculature.motility = 0.0035/60.0;//0.00035/60.0; //3.5/60.0
    std::cout<<"Vascular motility"<<vasculature.motility<<std::endl;
    vasculature.spatial_mechanical_factor=1.0;
    vasculature.spatial_proliferation_factor=1.0;
    vasculature.hypoxic_o2_threshold =1.0;
    vasculature.critical_o2_threshold = 1.0;
    vasculature.max_cells = 1.0;
    
    vasculature.other_properties_names.push_back("vascular_death_rate");
    vasculature.other_properties.push_back( 1.0/18/60 );  // 0
    vasculature.other_properties_names.push_back("vascular_proliferation_threshold");
    vasculature.other_properties.push_back( 0.001 );  // 1
    vasculature.other_properties_names.push_back("vascular_proliferation_saturation");
    vasculature.other_properties.push_back( 0.5 );   // 2
    vasculature.other_properties_names.push_back("vascular_chemotaxis_threshold");
    vasculature.other_properties.push_back( 0.001 );   // 3
    vasculature.other_properties_names.push_back("vascular_chemotaxis_saturation");
    vasculature.other_properties.push_back( 0.5 );   // 4
    vasculature.other_properties_names.push_back("effective_vascular_cutoff_threshold");
    vasculature.other_properties.push_back( 1.0E-8 );  // 5
    
    // These vectors store the base rates.  Environmental modifications made on the fly, just like the birth and death rates.  There is one three member vector (live, apop, necr) per substrate
    vasculature.secretion_rate[0][0]=9.9; // 1/time  Oxygen max secrection from vasculature
    vasculature.substrate_target[0][0]=1.0; // stuff (at times perhaps belongs in substrate/properties??

    
    // Initialize Population Structures with phenotype and physical (mesh) information
    
    for(int i=0; i<tissue.voxel_population_vectors.size(); i++)
     {
         
         // Tumor cells
         
         tissue.voxel_population_vectors[i].add_population("Johns_cancer", breast_cancer);
         
         // Makes the max population per voxle proportional to voxel volume (as required by voxels that are not of equal volume, like shells)
         tissue.voxel_population_vectors[i].phenotypes_vector[0].max_cells = tissue.tissue_mesh.voxels[i].volume/cell_volume;
         std::cout<<tissue.tissue_mesh.voxels[i].volume<<std::endl;
         std::cout<<tissue.tissue_mesh.voxels[i].volume/cell_volume<<std::endl;
         // Vasculature
         
         tissue.voxel_population_vectors[i].add_population("Vasculature", vasculature);
         
     }
    
    // Circle of vasculature and center voxel of "tumor"
    
//    tissue.voxel_population_vectors[84].live_cell_counts[0] = 16.0;
//
//    for(int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
//    {
//        if((tissue.tissue_mesh.voxels[i].center[0]*tissue.tissue_mesh.voxels[i].center[0]
//            +tissue.tissue_mesh.voxels[i].center[1]*tissue.tissue_mesh.voxels[i].center[1])>(tumor_radius-36)*(tumor_radius-36))
//        {
//            tissue.voxel_population_vectors[i].live_cell_counts[1] = 1.0;
//        }
//    }
//    
    // End circle of vasculature initialization
    
    // Parallel and opposing lines cells and vasculature - works for circular domain with radius 270 um

//    tissue.voxel_population_vectors[1].live_cell_counts[0] = 29.0;
//    tissue.voxel_population_vectors[2].live_cell_counts[0] = 29.0;
//    tissue.voxel_population_vectors[3].live_cell_counts[0] = 29.0;
//    tissue.voxel_population_vectors[4].live_cell_counts[0] = 29.0;
//    tissue.voxel_population_vectors[5].live_cell_counts[0] = 29.0;
//    tissue.voxel_population_vectors[6].live_cell_counts[0] = 29.0;
//    tissue.voxel_population_vectors[0].live_cell_counts[0] = 29.0;
////    tissue.voxel_population_vectors[151].live_cell_counts[1] = 1.0;
//    tissue.voxel_population_vectors[160].live_cell_counts[1] = 1.0;
//    tissue.voxel_population_vectors[159].live_cell_counts[1] = 1.0;
//    tissue.voxel_population_vectors[158].live_cell_counts[1] = 1.0;
//    tissue.voxel_population_vectors[162].live_cell_counts[1] = 1.0;
//    tissue.voxel_population_vectors[161].live_cell_counts[1] = 1.0;

    // end of opposing-line spacing
    
    // Spherical Initialization Pattern = live cells in center of domain and vasculature in outer voxel
    
    tissue.voxel_population_vectors[0].live_cell_counts[0] = 29.0;
    tissue.voxel_population_vectors[tissue.voxel_population_vectors.size()-1].live_cell_counts[1] = 0.5;
    
    // End spherical initialization
    
    int number_of_substrates = 2;

    tissue.create_substrate_vectors();
    
    for(int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
    {
        tissue.substrate_vectors[i].add_substrate("Oxygen", "mmHg");
        tissue.substrate_vectors[i].add_substrate("AF", "dimensionless");
        tissue.substrate_vectors[i].set_decay_constant ( 0, (0.01 * breast_cancer.uptake_rate[0][0]));
        tissue.substrate_vectors[i].set_decay_constant ( 1, 0.01667);  // 100-200 um
    }
    
    tissue.create_substrate_properties_vectors();  // Overall, these two components will ONLY handle diffusion.  Secretetion/decay/uptake are all separate.  They will need to be within the voxels - this stuff will ONLY be duffion, once the other stuff is set in a separate function.  This will be core stuff and the S/S will be specific (as it is in BioFVM. )
    
    Substrate_Properties substrate_properties;
    
    for(int i=0; i<tissue.tissue_mesh.voxel_links.size(); i++)
    {
        tissue.substrate_properties_vectors[i].add_diffusion_coefficient(10000);
    }
    
    tissue.tissue_mesh.write_mesh_to_matlab("./Output/mesh.mat");
    tissue.tissue_mesh.write_links_to_matlab("./Output/links.mat");
    tissue.write_voxel_populations_to_matlab("./Output/populations.mat");
    tissue.write_substrates_to_matlab("./Output/substrates.mat");
    tissue.write_substrate_properties_to_matlab("./Output/properties.mat");
    
    // Time, output, and simulation loop set up
   
    double t     = 0.0;
    double t_max = 1440.0*20;
    double dt    = 0.0025;//0.0025;
    std::cout<<tissue.formatted_minutes_to_DDHHMM( t )<<std::endl;
    double output_interval  = 60.0;  // how often you save data
    double next_output_time = 0.0;     // next time you save data
    double output_factor = output_interval/100; // Adjusts file output names for use with decimal output intervals (if file name would have a decimal in it, the file will be overwritten w/o this

    while (t < t_max)
     {
         // if it's time, save the simulation
         
         if( fabs( t - next_output_time ) < dt/2.0 )
         {
//             std::ofstream fileoutput4;
             char* filename;
             filename = new char [1024];

//             std::cout<<output_factor*next_output_time<<std::endl;
             sprintf( filename, "./Output/output_%06.0f.mat" , next_output_time*output_factor );
             
             tissue.write_all_to_matlab(filename, t);  //  Tumor - live, ap, nec; Vas - live, ap, nec; Oxygen; AF
             
             double height = 1500;
             double width = 1500;
             
             sprintf( filename, "./SVG/live_cells_%06.0f.svg" , next_output_time*output_factor );

             tissue.write_population_svg(0, t, filename, height, width);  // Live Tumor
//             tissue.write_tissue_apoptotic_svg(filename, height, width);  // apopotic tumor
//             tissue.write_tissue_necrotic_svg(filename, height, width);  // necrotic tumor
             
             sprintf( filename, "./SVG/vasculature_%06.0f.svg" , next_output_time*output_factor );
             tissue.write_population_svg(1, t, filename, height, width);  // Vasculature - needs moved to its own thing
             
             sprintf( filename, "./SVG/oxygen_%06.0f.svg" , next_output_time*output_factor );
             tissue.write_subsrate_svg(0, t, filename, height, width);  // Oxygen

             sprintf( filename, "./SVG/angiogenic_factor_%06.0f.svg" , next_output_time*output_factor );
             tissue.write_subsrate_svg(1, t, filename, height, width);  // AF
             
             next_output_time += output_interval;
         }

      for (int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
       {
           double temp;
           tissue.voxel_population_vectors[i].update_populations( dt, tissue.substrate_vectors[i].substrate_quantity[0]  );  // Oxygen dependent growth for cell populations.
           tissue.voxel_population_vectors[i].update_vascular_population ( dt, tissue.substrate_vectors[i].substrate_quantity[1] );
           temp = tissue.voxel_population_vectors[i].update_AF_secretion_rate( tissue.substrate_vectors[i].substrate_quantity[0] );  // o2 in voxel i
           tissue.voxel_population_vectors[i].phenotypes_vector[0].secretion_rate[1][0] = temp;
       }
         
         tissue.update_substrates( dt );
         tissue.flux_vascular_density( dt );
         tissue.flux_cell_populations( dt );
         tissue.run_diffusion ( dt, 0);
         tissue.run_diffusion ( dt, 1);
         t += dt;
         
     }

    // end output function
    
    RUNTIME_TOC();
    std::cout << std::endl << "Program wall time: ";
    display_stopwatch_value( std::cout , runtime_stopwatch_value() );
    std::cout << std::endl << "Done!" << std::endl;
    
    return 0;

}
    
