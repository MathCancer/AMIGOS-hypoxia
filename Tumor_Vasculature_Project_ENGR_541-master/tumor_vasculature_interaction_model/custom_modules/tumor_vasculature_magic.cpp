/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.3.0) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.3.0) [1],    #
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

#include "./tumor_vasculature_magic.h"
#include "./vasculature.h"
// declare cell definitions here 

Cell_Definition motile_cell; 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( time(NULL) ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = flow_cytometry_separated_cycle_model; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;
    cell_defaults.functions.custom_cell_rule = VEGF_secretion_and_vascular_death_function;
	
	// needed for a 2-D simulation: 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );
    int VEGF_substrate_index = microenvironment.find_density_index("VEGF");

	int G0G1_index = Ki67_advanced.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = Ki67_advanced.find_phase_index( PhysiCell_constants::S_phase );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

    // set default motiltiy
    cell_defaults.phenotype.motility.is_motile = false;
    cell_defaults.functions.update_migration_bias = oxygen_taxis_motility; // probably want to turn this off - don't want cells to move after all.
    cell_defaults.phenotype.motility.persistence_time = 15.0;
    cell_defaults.phenotype.motility.migration_bias = 0.5;
    cell_defaults.phenotype.motility.migration_speed = 0.05; // 0.5
    
	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10.0;
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38.0;
	
    // VEGF
    cell_defaults.phenotype.secretion.secretion_rates[VEGF_substrate_index] = 0.0; // DEFAULT STATE FOR A NORMOXIC CELL IS THAT THEY DON'T Secrete
    cell_defaults.phenotype.secretion.uptake_rates[VEGF_substrate_index] = 0.0;
    cell_defaults.phenotype.secretion.saturation_densities[VEGF_substrate_index] = 1.0;
    
    
    
    // add custom data here, if any
	
//    cell_defaults.custom_data.add_variable( "receptor" , "dimensionless", 0.0 );
//    cargo_cell.custom_data["receptor"] = 1.0;
//    pTemp->custom_data[ "receptor" ] = 0.0;
//    nearby[i]->custom_data["receptor"] = 0.0;

    
	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
//    motile_cell = cell_defaults;
//    motile_cell.type = 1;
//    motile_cell.name = "motile tumor cell";
//
//    // make sure the new cell type has its own reference phenotyhpe
//
//    motile_cell.parameters.pReference_live_phenotype = &( motile_cell.phenotype );
//
//    // enable random motility
//    motile_cell.phenotype.motility.is_motile = true;
//    motile_cell.phenotype.motility.persistence_time = 15.0; // 15 minutes
//    motile_cell.phenotype.motility.migration_speed = 0.25; // 0.25 micron/minute
//    motile_cell.phenotype.motility.migration_bias = 0.0;// completely random
//
//    // Set cell-cell adhesion to 5% of other cells
//    motile_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0.05;
//
//    // Set apoptosis to zero
//    motile_cell.phenotype.death.rates[apoptosis_model_index] = 0.0;
//
//    // Set proliferation to 10% of other cells.
//    // Alter the transition rate from G0G1 state to S state
//    motile_cell.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 0.1;
//
	return;
}

void setup_microenvironment( void )
{
    int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );
    
	// set domain parameters 
	
	default_microenvironment_options.X_range = {-500, 500}; 
	default_microenvironment_options.Y_range = {-500, 500}; 
	default_microenvironment_options.simulate_2D = true; // 2D! 
	
	// gradients need for this example (I think ...)

	default_microenvironment_options.calculate_gradients = false;
    
    // let BioFVM use oxygen as the default
    
    default_microenvironment_options.use_oxygen_as_first_field = true;
	
	// set Dirichlet conditions 

    default_microenvironment_options.outer_Dirichlet_conditions = false; // Do we need this?

    // SEE VASCULATURE.CPP --> Coarse Vasculature Setup = coarse_vasculature_setup
    //    default_microenvironment_options.Dirichlet_condition_vector[oxygen_substrate_index] = 10.0;  // Do we need this?
	
	// if there are more substrates, resize accordingly 
//    std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
//    default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
    
//    microenvironment.add_density( "VEGF", "dimensionless" , 3.5e3 , 0.09 );
	
    coarse_vasculature_setup();
    
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
    std::cout << "need to set up source functions!!!" << std::endl;
    
    for(int i = 0; i<coarse_vasculature.vascular_densities.size();i++)
    {
        std::cout<<coarse_vasculature.vascular_densities[i].functional<<std::endl;
    }
    
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;

//    pC = create_cell();
//    pC->assign_position( 0.0, 0.0, 0.0 );
//
//    pC = create_cell();
//    pC->assign_position( -100, 0, 0.0 );
//
//    pC = create_cell();
//    pC->assign_position( 0, 100, 0.0 );
//
//    // now create a motile cell
//
//    pC = create_cell(  );
//    pC->assign_position( 15.0, -18.0, 0.0 );

double cell_radius = cell_defaults.phenotype.geometry.radius;
double cell_spacing = 0.95 * 2.0 * cell_radius;

double tumor_radius = 125.0;

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
//            coarse_vasculature.vascular_densities
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

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with the Ki67 coloring 
	
	std::vector<std::string> output = false_cell_coloring_Ki67(pCell); 
	
	// if the cell is motile and not dead, paint it black 
	
	if( pCell->phenotype.death.dead == false && 
		pCell->type == 1 )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 	
	}
	
	return output; 
}

void VEGF_secretion_and_vascular_death_function(Cell* pCell, Phenotype& phenotype, double dt )
{
    
    // VEGF SECRETION CODE
    
    static int oxygen_index = pCell->get_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_substrate_index = pCell->get_microenvironment()->find_density_index("VEGF");
    
//    std::vector<double> substrates = pCell->nearest_density_vector();
//    double O2 = substrates[O2_i];
    double oxygen = (pCell->nearest_density_vector())[oxygen_index]; // does same thing as above two lines I think ...
    
    double r_base_VEGF = 10.0; // min-1
    double hypoxic_o2_threshold = 20.0; // pO2 mmHg --> how do I get this into a cell? Do I need to? See above in cell decks.
    double critical_o2_threshold = 5.0; // pO2 mmHg
    
    if( oxygen > hypoxic_o2_threshold )
    {
        pCell->phenotype.secretion.secretion_rates[VEGF_substrate_index] = 0.0;
    }
    
    else if ( oxygen < critical_o2_threshold)
    {
        pCell->phenotype.secretion.secretion_rates[VEGF_substrate_index] = r_base_VEGF;
    }
    
    else
    {
        pCell->phenotype.secretion.secretion_rates[VEGF_substrate_index] = r_base_VEGF*((hypoxic_o2_threshold - oxygen)/(hypoxic_o2_threshold - critical_o2_threshold));
    }

    // END VEGF SECRETION CODE
    
    // VASCULAR DEATH CODE
    
    double vascular_degradation_rate_per_cell = 1e-6;
    
    coarse_vasculature( pCell ).functional = coarse_vasculature( pCell ).functional/(1+dt*vascular_degradation_rate_per_cell);
    
    // END VASCULAR DEATH CODE
    
//    V - (V-1) = -dt*lambda*V
//    V(1+dt*lambda) = (V-1)
//    V=(V-1)/(1+dt*lambda)
    
    
}

void oxygen_taxis_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" );
    
    phenotype.motility.migration_bias_direction = pCell->nearest_gradient( oxygen_i );
    normalize( &(phenotype.motility.migration_bias_direction) ) ;
    
    return;
}
