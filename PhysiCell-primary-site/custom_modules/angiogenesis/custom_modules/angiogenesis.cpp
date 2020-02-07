/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
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

#include "./angiogenesis.h"
#include "../modules/PhysiCell_settings.h"

Cell_Definition stalk_cell;
Cell_Definition tip_cell;

void create_stalk_cell_types( void )
{
	stalk_cell = cell_defaults; 
	
	stalk_cell.name = "stalk cell";
	stalk_cell.type = 1;
        stalk_cell.phenotype.geometry.radius = 2*cell_defaults.phenotype.geometry.radius;
        stalk_cell.phenotype.volume.total = pow(stalk_cell.phenotype.geometry.radius,3)*4.188790204786391;

	// turn off proliferation; 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 	stalk_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	//int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index(0); 
	
	// set default uptake and secretion 
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
        static int VEGF_ID = microenvironment.find_density_index( "VEGF" ); // 1
	
	// oxygen 
	stalk_cell.phenotype.secretion.secretion_rates[oxygen_ID] = 1.0; 
	stalk_cell.phenotype.secretion.uptake_rates[oxygen_ID] = 0.0; 	
        stalk_cell.phenotype.secretion.saturation_densities[oxygen_ID] = 1.0;
        // VEGF 
	stalk_cell.phenotype.secretion.secretion_rates[VEGF_ID] = 0.0; 
	stalk_cell.phenotype.secretion.uptake_rates[VEGF_ID] = 1.0;	
        stalk_cell.phenotype.secretion.saturation_densities[VEGF_ID] = 1.0;
	
        // set apoptosis to survive 10 days (on average) 
        int apoptosis_index = stalk_cell.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	stalk_cell.phenotype.death.rates[apoptosis_index] = 0;
                
	
	// set functions 
	
	stalk_cell.functions.update_phenotype = NULL; 
        stalk_cell.functions.volume_update_function = NULL;
	stalk_cell.functions.custom_cell_rule = stalk_cell_rule; 
	stalk_cell.functions.update_migration_bias = NULL;	
	
        //set mechanics aspects
	tip_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 5.0;
	tip_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 5.0;
	
	return; 
}

void create_tip_cell_types( void )
{
	tip_cell = cell_defaults; 
	
	tip_cell.name = "tip cell";
	tip_cell.type = 2;
        tip_cell.phenotype.geometry.radius = 2*cell_defaults.phenotype.geometry.radius;
        tip_cell.phenotype.volume.total = pow(tip_cell.phenotype.geometry.radius,3)*4.188790204786391;

	// turn off proliferation; 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 	tip_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	//int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index(0); 
	
	// set default uptake and secretion 
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
        static int VEGF_ID = microenvironment.find_density_index( "VEGF" ); // 1
	
	// oxygen 
	tip_cell.phenotype.secretion.secretion_rates[oxygen_ID] = 1.0; 
	tip_cell.phenotype.secretion.uptake_rates[oxygen_ID] = 0.0; 	
        tip_cell.phenotype.secretion.saturation_densities[oxygen_ID] = 1.0;
        // VEGF 
	tip_cell.phenotype.secretion.secretion_rates[VEGF_ID] = 0.0; 
	tip_cell.phenotype.secretion.uptake_rates[VEGF_ID] = 1000.0;	
        tip_cell.phenotype.secretion.saturation_densities[VEGF_ID] = 1.0;
	
        // set apoptosis to survive 10 days (on average) 
        int apoptosis_index = tip_cell.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	tip_cell.phenotype.death.rates[apoptosis_index] = 0;
                
                
	// turn on motility; 
	tip_cell.phenotype.motility.is_motile = true; 
	tip_cell.phenotype.motility.persistence_time = 2; // 10.0; 
	tip_cell.phenotype.motility.migration_speed = 1; // 1;  
	tip_cell.phenotype.motility.migration_bias = 0.25; // 0.5;
	
	tip_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 5.0; // 0.0;
	tip_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 5.0; // 5.0;*/
	
	// set functions 
	
	tip_cell.functions.update_phenotype = NULL; 
        tip_cell.functions.volume_update_function = NULL;
	tip_cell.functions.custom_cell_rule = tip_cell_rule; 
	tip_cell.functions.update_migration_bias = tip_cell_motility;	
	
	// set custom data values 
	
	/*Parameter<double> paramD; 
	
	tip_cell.custom_data[ "oncoprotein" ] = 0.0; 
	tip_cell.custom_data[ "kill rate" ] = 
		parameters.doubles("tip_kill_rate"); // 1.0/15.0; // how often it tries to kill
	tip_cell.custom_data[ "attachment lifetime" ] = 
		parameters.doubles("tip_attachment_lifetime"); // 60.00; // how long it can stay attached 
	tip_cell.custom_data[ "attachment rate" ] = 
		parameters.doubles("tip_attachment_rate");*/ // 1.0/5.0; // how long it wants to wander before attaching	
	
	return; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints( "random_seed" ) ); 
	
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
	
	cell_defaults.parameters.o2_proliferation_saturation = 1.0;  
	cell_defaults.parameters.o2_reference = 1.0; 
	
	// set default uptake and secretion 
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
        static int VEGF_ID = microenvironment.find_density_index( "VEGF" ); // 1
	
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_ID] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_ID] = 0.1; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_ID] = 1.0;
        // VEGF 
	cell_defaults.phenotype.secretion.secretion_rates[VEGF_ID] = 1.0; 
	cell_defaults.phenotype.secretion.uptake_rates[VEGF_ID] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[VEGF_ID] = 1.0;

	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein; 
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// add custom data 
	
	cell_defaults.custom_data.add_variable( "oncoprotein" , "dimensionless", 1.0 ); 
        
        // create the endothelial cell type 
        create_stalk_cell_types();
        create_tip_cell_types();
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters

/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	// make sure ot override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/*
	All this is now in XML as of 1.6.0 
	
	// no gradients needed for this example 
	
	default_microenvironment_options.calculate_gradients = false; 
	
	// let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = true; 
	
	// set Dirichlet conditions 
	
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // normoxic conditions 
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/	
			
	initialize_microenvironment();


	return; 
}	

void introduce_stalk_cells( void )
{	
	// now seed stalk cells 
	
	int number_of_stalk_cells = 98; 
        double dx = 20.0, pos=-995.0;
	
	for( int i=0 ;i < number_of_stalk_cells ; i++ )
	{
		
		Cell* pCell = create_cell( stalk_cell ); 
		pCell->assign_position( -990, pos+dx, 0 ); 
                pos = pos + dx;
	}
	
	return; 
}

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	
	// Parameter<double> temp; 
	
	int i = parameters.doubles.find_index( "tumor_radius" ); 
	
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
	double p_mean = parameters.doubles( "oncoprotein_mean" ); 
	double p_sd = parameters.doubles( "oncoprotein_sd" ); 
	double p_min = parameters.doubles( "oncoprotein_min" ); 
	double p_max = parameters.doubles( "oncoprotein_max" ); 
	
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
			pCell->custom_data[0] = NormalRandom( p_mean, p_sd );
			if( pCell->custom_data[0] < p_min )
			{ pCell->custom_data[0] = p_min; }
			if( pCell->custom_data[0] > p_max )
			{ pCell->custom_data[0] = p_max; }
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );
				pCell->custom_data[0] = NormalRandom( p_mean, p_sd );
				if( pCell->custom_data[0] < p_min )
				{ pCell->custom_data[0] = p_min; }
				if( pCell->custom_data[0] > p_max )
				{ pCell->custom_data[0] = p_max; }				
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );
				pCell->custom_data[0] = NormalRandom( p_mean, p_sd );
				if( pCell->custom_data[0] < p_min )
				{ pCell->custom_data[0] = p_min; }
				if( pCell->custom_data[0] > p_max )
				{ pCell->custom_data[0] = p_max; }
		
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );
					pCell->custom_data[0] = NormalRandom( p_mean, p_sd );
					if( pCell->custom_data[0] < p_min )
					{ pCell->custom_data[0] = p_min; }
					if( pCell->custom_data[0] > p_max )
					{ pCell->custom_data[0] = p_max; }
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	double sum = 0.0; 
	double min = 9e9; 
	double max = -9e9; 
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data[0]; 
		sum += r;
		if( r < min )
		{ min = r; } 
		if( r > max )
		{ max = r; }
	}
	double mean = sum / ( all_cells->size() + 1e-15 ); 
	// compute standard deviation 
	sum = 0.0; 
	for( int i=0; i < all_cells->size(); i++ )
	{
		sum +=  ( (*all_cells)[i]->custom_data[0] - mean )*( (*all_cells)[i]->custom_data[0] - mean ); 
	}
	double standard_deviation = sqrt( sum / ( all_cells->size() - 1.0 + 1e-15 ) ); 
        
	
	std::cout << std::endl << "Oncoprotein summary: " << std::endl
			  << "===================" << std::endl; 
	std::cout << "mean: " << mean << std::endl; 
	std::cout << "standard deviation: " << standard_deviation << std::endl; 
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl; 
	
	return; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" );
    
        update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
        //if (pCell->type == 1)
	//std::cout << "Cell: "<< pCell->type_name << "Type: " << pCell->type << std::endl;
	// if cell is dead, don't bother with future phenotype changes. 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}
        
        // multiply proliferation rate by the oncoprotein 
	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}

std::vector<std::string> angiogenesis_coloring_function( Cell* pCell )
{
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	
	static double p_min = parameters.doubles( "oncoprotein_min" ); 
	static double p_max = parameters.doubles( "oncoprotein_max" ); 
	
	// immune are black
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
            output[0] = "red";
            output[1] = "black";
            output[2] = "red"; 
            return output; 
        }
        
	
	// live cells are green, but shaded by oncoprotein value 
	if( pCell->phenotype.death.dead == false )
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->custom_data[oncoprotein_i]-p_min) * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
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

void tip_cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if attached, biased motility towards director chemoattractant 
	// otherwise, biased motility towards cargo chemoattractant 
	
	static int endothelial_factor_index = microenvironment.find_density_index( "VEGF" ); 

	// if not docked, attempt biased chemotaxis 
	if( pCell->state.neighbors.size() == 0 )
	{
		phenotype.motility.is_motile = false; 			
	}
	else
	{
		phenotype.motility.is_motile = true; 
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(endothelial_factor_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );
	}
	
	return; 
}

void tip_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    for(int i =0; i < pCell->state.neighbors.size();i++ )
    {
        //transition tip to stalk
        if ((pCell->state.neighbors[i])->type_name == "tip cell")
        {
            pCell->convert_to_cell_definition(tip_cell);
        }
    }
    
}

void stalk_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int endothelial_factor_index = microenvironment.find_density_index( "VEGF" );
  
    if (pCell->nearest_density_vector()[endothelial_factor_index] > 0.12)
    {
        printf("VEGF concentration: %e -- Cell: ",pCell->nearest_density_vector()[endothelial_factor_index]);
        std::cout << pCell->type_name << std::endl;
        Cell* NewCell = create_cell( tip_cell ); 
        NewCell->assign_position( pCell->position[0]+2*pCell->phenotype.geometry.radius, pCell->position[1], 0 ); 
    }
    
}
