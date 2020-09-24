/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project             #
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

#include "./hypoxia2D.h"

Cell_Definition stalk_cell;
Cell_Definition blood_cell;

void create_stalk_cell_types( void )
{
	stalk_cell = cell_defaults; 
	
	stalk_cell.name = "stalk cell";
	stalk_cell.type = 1;
        stalk_cell.phenotype.geometry.radius = 2*cell_defaults.phenotype.geometry.radius;
        stalk_cell.phenotype.volume.total = pow(stalk_cell.phenotype.geometry.radius,3)*4.188790204786391;

	// turn off proliferation; 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 	

	stalk_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	//int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index(0); 
	
	// set default uptake and secretion 
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
        static int VEGF_ID = microenvironment.find_density_index( "VEGF" ); // 1
	
	// oxygen 
	stalk_cell.phenotype.secretion.secretion_rates[oxygen_ID] = 90.0; 
	stalk_cell.phenotype.secretion.uptake_rates[oxygen_ID] = 0.0; 	
        stalk_cell.phenotype.secretion.saturation_densities[oxygen_ID] = 90.0;
        // VEGF 
	stalk_cell.phenotype.secretion.secretion_rates[VEGF_ID] = 0.0; 
	stalk_cell.phenotype.secretion.uptake_rates[VEGF_ID] = 1.0;	
        stalk_cell.phenotype.secretion.saturation_densities[VEGF_ID] = 1.0;
	
        // turn off apoptosis
        int apoptosis_index = stalk_cell.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	stalk_cell.phenotype.death.rates[apoptosis_index] = 0;
                
	
	// set functions 
	
	stalk_cell.functions.update_phenotype = stalk_cell_rule; 
        stalk_cell.functions.volume_update_function = NULL;
	//stalk_cell.functions.custom_cell_rule = stalk_cell_rule; 
	stalk_cell.functions.custom_cell_rule = NULL; 
	stalk_cell.functions.update_migration_bias = NULL;
	stalk_cell.functions.update_velocity = NULL;	
	
        //set mechanics aspects
	stalk_cell.phenotype.motility.is_motile = false; 
	stalk_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0.0;
	stalk_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 100.0;

	//stalk_cell.functions.custom_cell_rule = extra_elastic_attachment_mechanics; 

	// tip cells ID
	//stalk_cell.custom_data.add_variable( "tip_cell ID" , 0.0); 
	
	return; 
}

void create_blood_cell_types( void )
{
	blood_cell = cell_defaults; 
	
	blood_cell.name = "blood cell";
	blood_cell.type = 2;
        blood_cell.phenotype.geometry.radius = cell_defaults.phenotype.geometry.radius;
        blood_cell.phenotype.volume.total = pow(stalk_cell.phenotype.geometry.radius,3)*4.188790204786391;

	// turn off proliferation; 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 	

	blood_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	//int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index(0); 
	
	// set default uptake and secretion 
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
        static int VEGF_ID = microenvironment.find_density_index( "VEGF" ); // 1
	
	// oxygen 
	blood_cell.phenotype.secretion.secretion_rates[oxygen_ID] = 90.0; 
	blood_cell.phenotype.secretion.uptake_rates[oxygen_ID] = 0.0; 	
    blood_cell.phenotype.secretion.saturation_densities[oxygen_ID] = 90.0;
    // VEGF 
	blood_cell.phenotype.secretion.secretion_rates[VEGF_ID] = 0.0; 
	blood_cell.phenotype.secretion.uptake_rates[VEGF_ID] = 0.0;	
    blood_cell.phenotype.secretion.saturation_densities[VEGF_ID] = 0.0;
	
    // turn off apoptosis
    int apoptosis_index = blood_cell.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	blood_cell.phenotype.death.rates[apoptosis_index] = 0;
                
	
	// set functions 
	
	blood_cell.functions.update_phenotype = NULL; 
        blood_cell.functions.volume_update_function = NULL;
	blood_cell.functions.custom_cell_rule = NULL; 
	blood_cell.functions.update_migration_bias = NULL;	
	
        //set mechanics aspects
	blood_cell.phenotype.motility.is_motile = false; 
	blood_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 1.0;
	blood_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 1.0;

	blood_cell.functions.custom_cell_rule = NULL; 

	// tip cells ID
	//stalk_cell.custom_data.add_variable( "tip_cell ID" , 0.0); 
	
	return; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom(0); 
	
	int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" ); 
	int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" ); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// use default proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.parameters.o2_proliferation_saturation =  38.0;  
	cell_defaults.parameters.o2_reference = cell_defaults.parameters.o2_proliferation_saturation; 
	
	// set default motiltiy
	cell_defaults.phenotype.motility.is_motile = false; 
	cell_defaults.functions.update_migration_bias = oxygen_taxis_motility; 
	cell_defaults.phenotype.motility.persistence_time = 15.0; 
	cell_defaults.phenotype.motility.migration_bias = 0.5; 
	cell_defaults.phenotype.motility.migration_speed = 0.05; // 0.5 
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_i] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_i] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_i] = 38; 
	
	// VEGF 
	cell_defaults.phenotype.secretion.secretion_rates[VEGF_i] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[VEGF_i] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[VEGF_i] = 1.0; 

	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype; 
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// add custom data 
	
	std::vector<double> genes = { 1.0, 0.0 }; // RFP, GFP 
	std::vector<double> proteins = {1.0, 0.0 }; // RFP, GFP; 
	
	double default_degradation_rate = 6.8e-5; // 7-day half-life 
	// 4.8e-4; // 24 hour half-life 
	// 0.0077; // 90 minute half-life 
	// 0.019; // 90% degrades in 120 minutes 
	
	std::vector<double> degradation_rates = { default_degradation_rate , default_degradation_rate }; 
	
	double default_production_rate = 4.8e-4; // 24 hour half ramp-up 
	// 0.0019; // 6 hours to reach 50% 
	// 0.0068; // 1.7 hours to reach 50% 
	// 0.23; // 10 minutes to reach 90% 
	
	std::vector<double> creation_rates = { default_production_rate , default_production_rate }; 
	
	cell_defaults.custom_data.add_vector_variable( "genes" , "dimensionless", genes ); 
	
	cell_defaults.custom_data.add_vector_variable( "proteins" , "dimensionless", proteins ); 
	cell_defaults.custom_data.add_vector_variable( "creation_rates" , "1/min" , creation_rates ); 
	cell_defaults.custom_data.add_vector_variable( "degradation_rates" , "1/min" , degradation_rates ); 

	cell_defaults.custom_data.add_variable( "persistence time" , "dimensionless" , 0.0 ); 
	
	cell_defaults.custom_data.add_variable( "hypoxic memory" , "min" , 0.0 ); 
	
	
	std::vector<double> color = {255, 255, 255};
	cell_defaults.custom_data.add_vector_variable( "nuclear_color" , "dimensionless", color ); 
	cell_defaults.custom_data.add_vector_variable( "cytoplasmic_color" , "dimensionless", color ); 

	//Mechanics
	Parameter<double> paramD; 
	
	// for cargo-worker 
	paramD = parameters.doubles["elastic_coefficient"]; 
	cell_defaults.custom_data.add_variable( "elastic coefficient" , paramD.units, paramD.value ); 

	// create the endothelial cell type 
        create_stalk_cell_types();
	create_blood_cell_types();
	
	return; 
}

void setup_microenvironment( void )
{
	default_microenvironment_options.X_range = {-1500, 1500}; 
	default_microenvironment_options.Y_range = {-1500, 1500}; 
	// make sure ot override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}

	initialize_microenvironment(); 	

	return; 
}	

void introduce_stalk_cells( void )
{	
	// now seed stalk cells 
	
	int number_of_stalk_cells = 148; 
        double dx = 20.0, pos=-1475.0;

	double radius = 100;
	double theta = 0;
	const double PI = 3.14159265;
	for( int i=0 ;i < 40; i++ )
	{	
		Cell* pCell = create_cell( stalk_cell ); 
		pCell->assign_position( radius*cos(theta) + 750, radius*sin(theta)+ 750, 0 ); 
                theta = theta + 9*(PI/180);
	}

	/*for( int i=0 ;i < number_of_stalk_cells ; i++ )
	{
		
		Cell* pCell = create_cell( stalk_cell ); 
		pCell->assign_position( -1495.0, pos, 0 ); 
                pos = pos + dx;
	}*/
	/*for( int i=0 ;i < number_of_stalk_cells ; i++ )
	{
		
		Cell* pCell = create_cell( stalk_cell ); 
		pCell->assign_position( 1495.0, pos, 0 ); 
                pos = pos - dx;
	}*/
	
	//Neighboors
	for(int i=0;i<all_cells->size();i++)
	{
	  for(int j=i;j<all_cells->size();j++)
	  {
	    if ((*all_cells)[i]->type == 0 || (*all_cells)[j]->type == 0 || i == j) continue;
	    if (dist((*all_cells)[i]->position,(*all_cells)[j]->position) < (*all_cells)[i]->phenotype.geometry.radius)
	    {
	      (*all_cells)[i]->state.neighbors.push_back((*all_cells)[j]);
	      (*all_cells)[j]->state.neighbors.push_back((*all_cells)[i]);
	      printf("ID: %d - Neighbors: %d -- tipID of the first neighbor: %d\n",i,(*all_cells)[i]->state.neighbors.size(),j);
	    }	
	  }
  	}

	double cell_radius = blood_cell.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double vessel_radius = radius-10; 
		
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = vessel_radius; 
	double y = 0.0; 
	
	int n = 0; 
	while( y < vessel_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( vessel_radius*vessel_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(blood_cell); // tumor cell 
			pCell->assign_position( 750.0+x , 750.0+y , 0.0 );
					
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(blood_cell); // tumor cell 
				pCell->assign_position( 750.0+x , 750.0-y , 0.0 );
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(blood_cell); // tumor cell 
				pCell->assign_position( 750.0-x , 750.0+y , 0.0 );
								
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(blood_cell); // tumor cell 
					pCell->assign_position( 750.0-x , 750.0-y , 0.0 );
					
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	return; 
}

void setup_tissue( void )
{
	static int genes_i = 0; 
	static int proteins_i =1; 
	static int creation_rates_i = 2; 
	static int degradation_rates_i = 3; 
	
	static int red_i = 0; 
	static int green_i = 1; 	
	
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = 250.0; 
		
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
	int n = 0; 
	while( y < tumor_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(); // tumor cell 
			pCell->assign_position( x , y , 0.0 );
					
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );
								
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );
					
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	//To old models comment this
	//introduce_stalk_cells();
		
	return; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	//std::cout << "Position: ( " << pCell->position[0] << " " << pCell->position[1] <<" " << pCell->position[2] << ") -- (" << "type: " << pCell->type << std::endl;
	static int genes_i = 0; 
	static int proteins_i =1; 
	static int creation_rates_i = 2; 
	static int degradation_rates_i = 3; 
	
	static int red_i = 0; 
	static int green_i = 1; 	

	static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" ); 
	static int persistence_time_i = pCell->custom_data.find_variable_index( "persistence time" );
	
	static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" ); 
	static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" ); 
	
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 
	
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
//	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 

	// set genes 
	
	double pO2 = (pCell->nearest_density_vector())[oxygen_i]; 
	
	static double FP_hypoxic_switch = 10.0; // 
	static double phenotype_hypoxic_switch = 10.0;  // 
	
	// permanent gene switch 
	if( pO2 < FP_hypoxic_switch )
	{
		pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0; 
		pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0; 
	}
	
	if( pO2 < phenotype_hypoxic_switch )
	{
		phenotype.secretion.secretion_rates[VEGF_i] = 10.0; 
	}
	else
	{
		phenotype.secretion.secretion_rates[VEGF_i] = 0.0; 
	}

	// update the proteins
	
	for( int i=0; i < pCell->custom_data.vector_variables[genes_i].value.size(); i++ )
	{
		double temp = pCell->custom_data.vector_variables[creation_rates_i].value[i]; // alpha_i
		temp += pCell->custom_data.vector_variables[degradation_rates_i].value[i]; // alpha_i + beta_i 
		temp *= pCell->custom_data.vector_variables[genes_i].value[i]; // G_i^n ( alpha_i + beta_i ); 
		temp *= dt; // dt*G_i^n ( alpha_i + beta_i ); 
		pCell->custom_data.vector_variables[proteins_i].value[i] += temp; // P_i = P_i + dt*G_i^n ( alpha_i + beta_i ); 
		temp = pCell->custom_data.vector_variables[creation_rates_i].value[i]; // alpha_i 
		temp *= pCell->custom_data.vector_variables[genes_i].value[i]; // G_i^n * alpha_i 
		temp += pCell->custom_data.vector_variables[degradation_rates_i].value[i]; // G_i^n * alpha_i + beta_i 
		temp *= dt; // dt*( G_i^n * alpha_i + beta_i ); 
		temp += 1.0; // 1.0 + dt*( G_i^n * alpha_i + beta_i ); 
		pCell->custom_data.vector_variables[proteins_i].value[i] /= temp; // P_i = ( P_i + dt*G_i^n ( alpha_i + beta_i ) ) / ( 1.0 + dt*( G_i^n * alpha_i + beta_i ) ); 
	}
	

	
	if( pO2 < phenotype_hypoxic_switch && pCell->type != 1 )
	{
		if (phenotype.motility.is_motile == false) pCell->custom_data[persistence_time_i] = 0.0;
		phenotype.motility.is_motile = true; 
		phenotype.motility.migration_speed = 0.25; //   migration_bias = 0.85; 
		// phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1; 
	}

	// Model 1 -- No phenotypic persistence
	/*else
	{
		phenotype.motility.migration_speed = 0.0;
	}*/

	// Model 2 -- Phenotypic permanence
	// No change

	// Model 3 -- Phenotypic persistence
	/*else
	{
		static double persistence_time = 1440; // // 5760; // 4 days // 3600; // 60 hours // 360.0; // 6 hours // too short! 
		static double probability = dt/persistence_time; 
		
		if( phenotype.motility.is_motile == true )
		{
			if( UniformRandom() < probability )
			{
				phenotype.motility.is_motile = false; 

			}
			else
			{
				// phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1; 
			}
			
		}
	}*/
	// Determinist version
	else
	{
		static double persistence_time = 6;
		pCell->custom_data[persistence_time_i]+= dt;	
		if (phenotype.motility.is_motile == true && pCell->custom_data[persistence_time_i] < persistence_time)
		{
			phenotype.motility.is_motile = false;
			pCell->custom_data[persistence_time_i] = 0.0;
		} 	
	}
	return; 
}

std::vector<std::string> AMIGOS_coloring_function( Cell* pCell )
{
	static int genes_i = 0; 
	static int proteins_i =1; 
	static int creation_rates_i = 2; 
	static int degradation_rates_i = 3; 
	
	static int red_i = 0; 
	static int green_i = 1; 

	// immune are black
	std::vector< std::string > output( 4, "black" ); 
/*	
	if( pCell->type == 1 )
	{ return output; } 
*/

	static int cyto_color_i = 4;
	static int nuclear_color_i = 5;

	
	// live cells are a combination of red and green 
	if( pCell->phenotype.death.dead == false )
	{
		int red   = (int) round( pCell->custom_data.vector_variables[proteins_i].value[red_i] * 255.0 ); 
		int green = (int) round( pCell->custom_data.vector_variables[proteins_i].value[green_i] * 255.0); 

		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,0)", red, green );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		
		pCell->custom_data.vector_variables[cyto_color_i].value[0] = red; 
		pCell->custom_data.vector_variables[cyto_color_i].value[1] = green; 
		pCell->custom_data.vector_variables[cyto_color_i].value[2] = 0.0; 
		
		pCell->custom_data.vector_variables[nuclear_color_i].value[0] = red / 2.0; 
		pCell->custom_data.vector_variables[nuclear_color_i].value[1] = green / 2.0; 
		pCell->custom_data.vector_variables[nuclear_color_i].value[2] = 0.0 / 2.0; 
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
		pCell->custom_data.vector_variables[cyto_color_i].value[0] = 255; 
		pCell->custom_data.vector_variables[cyto_color_i].value[1] = 0; 
		pCell->custom_data.vector_variables[cyto_color_i].value[2] = 0.0; 
		
		pCell->custom_data.vector_variables[nuclear_color_i].value[0] = 125; 
		pCell->custom_data.vector_variables[nuclear_color_i].value[1] = 0; 
		pCell->custom_data.vector_variables[nuclear_color_i].value[2] = 0; 		
		
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
		
		pCell->custom_data.vector_variables[cyto_color_i].value[0] = 250; 
		pCell->custom_data.vector_variables[cyto_color_i].value[1] = 138; 
		pCell->custom_data.vector_variables[cyto_color_i].value[2] = 38; 
		
		pCell->custom_data.vector_variables[nuclear_color_i].value[0] = 139; 
		pCell->custom_data.vector_variables[nuclear_color_i].value[1] = 69; 
		pCell->custom_data.vector_variables[nuclear_color_i].value[2] = 19;
	}	
	
	return output; 
}

void oxygen_taxis_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( oxygen_i ); 
	normalize( &(phenotype.motility.migration_bias_direction) ) ; 
	
	return; 
}

void stalk_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int endothelial_factor_index = microenvironment.find_density_index( "VEGF" );
    static int TipIDindex = 0;

    //std::cout << "Position: ( " << pCell->position[0] << " " << pCell->position[1] <<" " << pCell->position[2] << ") -- (" << "type: " << pCell->type << std::endl;
  
    /*if (pCell->nearest_density_vector()[endothelial_factor_index] > 0.12)
    {
        Cell* NewCell = create_cell( tip_cell );
        TipIDCount+= 1;
        NewCell->custom_data[TipIDindex] = TipIDCount;
        printf("VEGF concentration: %e -- ID: %e -- Total: %d\n",pCell->nearest_density_vector()[endothelial_factor_index],NewCell->custom_data[TipIDindex], all_cells->size());    
        NewCell->assign_position( pCell->position[0]+2*pCell->phenotype.geometry.radius, pCell->position[1], 0 );
    }*/
    
}

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	
	// dettach cells if too far apart 
	static double max_elastic_displacement = parameters.doubles("max_elastic_displacement");
	static double max_displacement_squared = max_elastic_displacement*max_elastic_displacement; 
	
	/*if( norm_squared( displacement ) > max_displacement_squared )
	{
		dettach_cells( pActingOn , pAttachedTo );
		std::cout << "\t\tDETACH!!!!!" << std::endl; 
		return; 
	}*/
	
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	//if (pCell->state.neighbors.size() != 0) printf("Neighbors: %d",pCell->state.neighbors.size());
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}
