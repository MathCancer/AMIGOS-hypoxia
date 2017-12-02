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

#include "./primary_site.h"

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom(0); 
	
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
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	// set default motiltiy
	cell_defaults.phenotype.motility.is_motile = false; 
	cell_defaults.functions.update_migration_bias = oxygen_taxis_motility; 
	cell_defaults.phenotype.motility.persistence_time = 15.0; 
	cell_defaults.phenotype.motility.migration_bias = 0.5; 
	cell_defaults.phenotype.motility.migration_speed = 0.5; 
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 38; 

	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype; 
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// add custom data 
	
	std::vector<double> genes = { 1.0, 0.0 }; // RFP, GFP 
	std::vector<double> proteins = {1.0, 0.0 }; // RFP, GFP; 
	
	double default_degradation_rate = 4.8e-4; // 24 hour half-life 
	// 0.0077; // 90 minute half-life 
	// 0.019; // 90% degrades in 120 minutes 
	
	std::vector<double> degradation_rates = { default_degradation_rate , default_degradation_rate }; // degrade by 90% in 120 minutes 
	
	double default_production_rate = 0.0019; // 6 hours to reach 50% 
	// 0.0068; // 1.7 hours to reach 50% 
	// 0.23; // 10 minutes to reach 90% 
	
	std::vector<double> creation_rates = { default_production_rate , default_production_rate }; // 10 minute time scale 
	
	cell_defaults.custom_data.add_vector_variable( "genes" , "dimensionless", genes ); 
	cell_defaults.custom_data.add_vector_variable( "proteins" , "dimensionless", proteins ); 
	cell_defaults.custom_data.add_vector_variable( "creation_rates" , "1/min" , creation_rates ); 
	cell_defaults.custom_data.add_vector_variable( "degradation_rates" , "1/min" , degradation_rates ); 
	
	cell_defaults.custom_data.add_variable( "hypoxic memory" , "min" , 0.0 ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters

	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
	
	// no gradients needed for this example 
	
	default_microenvironment_options.calculate_gradients = true; 
	
	// let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = true; 
	
	// set Dirichlet conditions 
	
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_condition_vector[0] = 60; // 38; // physioxic conditions 
		
	// add ECM 
	
	Microenvironment* pME = get_default_microenvironment(); 
	int ECM_i = pME->find_density_index( "ECM" ); 
	if( ECM_i < 0 )
	{
		std::cout << "Adding ECM to the microenvironment ... " << std::endl; 
		pME->add_density( "ECM", "dimensionless" , 0.0 , 0.0 ); 
		ECM_i = pME->find_density_index( "ECM" ); 
		
		default_microenvironment_options.Dirichlet_condition_vector[ECM_i] = 1;  
		default_microenvironment_options.Dirichlet_activation_vector[ECM_i] = false;
	}
	
	coarse_vasculature_setup(); // ANGIO 
		
	initialize_microenvironment(); 	

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
/*			
			pCell->custom_data[0] = NormalRandom( 1.0, 0.33 );
			if( pCell->custom_data[0] < 0.0 )
			{ pCell->custom_data[0] = 0.0; }
			if( pCell->custom_data[0] > 2.0 )
			{ pCell->custom_data[0] = .0; }
*/		
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );
/*				
				pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
				if( pCell->custom_data[0] < 0.0 )
				{ pCell->custom_data[0] = 0.0; }
				if( pCell->custom_data[0] > 2.0 )
				{ pCell->custom_data[0] = .0; }
*/			
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );
/*				
				pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
				if( pCell->custom_data[0] < 0.0 )
				{ pCell->custom_data[0] = 0.0; }
				if( pCell->custom_data[0] > 2.0 )
				{ pCell->custom_data[0] = .0; }
*/				
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );
/*					
					pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
					if( pCell->custom_data[0] < 0.0 )
					{ pCell->custom_data[0] = 0.0; }
					if( pCell->custom_data[0] > 2.0 )
					{ pCell->custom_data[0] = .0; }
*/				
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	return; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int genes_i = 0; 
	static int proteins_i =1; 
	static int creation_rates_i = 2; 
	static int degradation_rates_i = 3; 
	
	static int red_i = 0; 
	static int green_i = 1; 	

	static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" ); 
	
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
	
	static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	double pO2 = (pCell->nearest_density_vector())[oxygen_i]; 
	
	static double hypoxic_switch = 10.0; 
	
	if( pO2 < hypoxic_switch )
	{
		pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0; 
		pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0; 
	}

	// update the proteins
	
	for( int i=0; i < pCell->custom_data.vector_variables.size(); i++ )
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
	
	// set phenotype in response to temporary or permanent changes 
	
	// model 1
	// if hypoxic, motile. 
	if( pO2 < hypoxic_switch )
	{
		phenotype.motility.is_motile = true; 
	}
	else
	{
		phenotype.motility.is_motile = false; 
	}
	
	// model 2
	// if green, motile 
/*	
	if( pCell->custom_data.vector_variables[proteins_i].value[green_i] > 0.5 )
	{
		phenotype.motility.is_motile = true; 
	}
	else
	{
		phenotype.motility.is_motile = false; 
	}
*/
	
	// model 3
	// if green, motile. but only for awhile
/*	
	if( pO2 < hypoxic_switch )
	{ 
		pCell->custom_data[hypoxic_memory_i] += 3.0*dt; 
	}
	else
	{
		pCell->custom_data[hypoxic_memory_i] -= dt; 
		if( pCell->custom_data[hypoxic_memory_i] < 0.0 )
		{ 
			pCell->custom_data[hypoxic_memory_i] = 0.0; 
		}
	}
	if( pCell->custom_data.vector_variables[proteins_i].value[green_i] > 0.5 && pCell->custom_data[hypoxic_memory_i] > 0.0 )
	{
		phenotype.motility.is_motile = true; 
	}
	else
	{
		phenotype.motility.is_motile = false; 
	}
*/	
	
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
	
	if( pCell->type == 1 )
	{ return output; } 
	
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
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
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
