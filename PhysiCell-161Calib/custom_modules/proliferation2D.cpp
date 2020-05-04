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

#include "./proliferation2D.h"


void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom(0); 
	
	int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	/// turn off proliferation and apoptosis; 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	cell_defaults.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0;
	
	cell_defaults.parameters.o2_proliferation_saturation =  38.0;  
	cell_defaults.parameters.o2_reference = cell_defaults.parameters.o2_proliferation_saturation; 
	
	// set default motiltiy
	cell_defaults.phenotype.motility.is_motile = false;  
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_i] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_i] = 0.4663;
	
	/* cell_defaults.parameters.o2_proliferation_saturation = parameters.doubles["prol_saturation"].value; 
	cell_defaults.parameters.o2_proliferation_threshold = parameters.doubles["prol_threshold"].value; 
	
	cell_defaults.parameters.o2_reference = cell_defaults.parameters.o2_proliferation_saturation; */
	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype;
	
	
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
	cell_defaults.phenotype.mechanics.cell_BM_adhesion_strength = 0.0;
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
	cell_defaults.phenotype.mechanics.cell_BM_repulsion_strength = 0.0;
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0;
	
	cell_defaults.custom_data.add_variable( "PercentagemProl" , "dimensionless" , 0.0 );
	
	return; 
}

void setup_microenvironment( void )
{
	// make sure ot override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}	
	initialize_microenvironment();
	
	for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ )
		microenvironment.add_dirichlet_node(n,default_microenvironment_options.Dirichlet_condition_vector);

	return; 
}	

void setup_tissue( void )
{		
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = 4300.0; 
		
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
		
	return; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	//Update dirichlet nodes
	microenvironment.remove_dirichlet_node(microenvironment.nearest_voxel_index( pCell->position));
	
	/* double pO2 = (pCell->nearest_density_vector())[0]; 
	double multiplier = 1.0;
     if( pO2 < pCell->parameters.o2_proliferation_saturation )
     {
         multiplier = ( pO2 - pCell->parameters.o2_proliferation_threshold ) 
             / ( pCell->parameters.o2_proliferation_saturation - pCell->parameters.o2_proliferation_threshold );
     }
     if( pO2 < pCell->parameters.o2_proliferation_threshold )
     { 
         multiplier = 0.0; 
     }
     
	pCell->custom_data[0] = multiplier; */
	
	return; 
}

std::vector<std::string> AMIGOS_coloring_function( Cell* pCell )
{
	int red   = 255.0; 
	int green = 0; 
	std::vector< std::string > output( 4, "black" ); 
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,%u,0)", red, green );
	output[0].assign( szTempString );
	output[1].assign( szTempString );

	sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
	output[2].assign( szTempString );
		
	return output;
}

void oxygen_taxis_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( oxygen_i ); 
	normalize( &(phenotype.motility.migration_bias_direction) ) ; 
	
	return; 
}

void QOI(const char* FILE)
{
	std::ofstream OutFile (FILE);
	std::vector<double> KI67(43);
	
	std::vector< double > position = {0.0,0.0,0.0};
	/* std::vector< double > TempDistance(KI67.size(),100);
	std::vector< int > index(KI67.size(),0);
	
	for(int i=0;i<all_cells->size();i++){
		position[1] = 100;
		for(int j=0;j < KI67.size();j++){
			double distance = dist(position,(*all_cells)[i]->position);
			if(distance < TempDistance[j]){
				TempDistance[j] = distance;
				index[j] = i;
			}
			position[1] += 100;
		}
	} 
	for(int j=0;j < KI67.size();j++)
		OutFile <<  (j+1)*100 << "\t" << (*all_cells)[index[j]]->custom_data[0] <<std::endl;
	OutFile.close();*/

	for(int j=0;j < KI67.size();j++){
		position[1] += 100;
		std::vector< double > Value = microenvironment.nearest_density_vector( position );
		double pO2 = Value[0]; 
		double multiplier = 1.0;
		if( pO2 < parameters.doubles["prol_saturation"].value )
		{
			multiplier = ( pO2 - parameters.doubles["prol_threshold"].value ) 
				 / ( parameters.doubles["prol_saturation"].value - parameters.doubles["prol_threshold"].value );
		}
		if( pO2 < parameters.doubles["prol_threshold"].value )
		{ 
			multiplier = 0.0; 
		}
		 
		OutFile <<  (j+1)*0.1 << "\t" << multiplier <<  "\t" << pO2 << std::endl;
	}
	return; 
}
