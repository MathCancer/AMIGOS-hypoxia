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
#include "./vasculature.h"



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
    
    cell_defaults.phenotype.cycle.sync_to_cycle_model( Ki67_advanced );
    
    // Make sure we're ready for 2D
    
    cell_defaults.functions.set_orientation = up_orientation;
    cell_defaults.phenotype.geometry.polarity = 1.0;
    cell_defaults.phenotype.motility.restrict_to_2D = true;
    
    // use default proliferation and death
    
    int Ki67_negative_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_negative );
    int Ki67_positive_premitotic_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_positive_premitotic );
    int Ki67_positive_postmitotic_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_positive_postmitotic );
    
    cell_defaults.phenotype.cycle.data.transition_rate( Ki67_negative_index ,Ki67_positive_premitotic_index ) = parameters.doubles("rQ");  //0.004605;
    cell_defaults.phenotype.cycle.data.transition_rate( Ki67_positive_premitotic_index ,Ki67_positive_postmitotic_index ) = parameters.doubles("r1");  //0.001282;
    cell_defaults.phenotype.cycle.data.transition_rate( Ki67_positive_postmitotic_index, Ki67_negative_index ) = parameters.doubles("r2");  //0.0066667;
    
    int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
    
    cell_defaults.phenotype.death.rates[apoptosis_index] = parameters.doubles( "rA" ); //0.0000387
    
    
    
    cell_defaults.parameters.o2_proliferation_saturation =  parameters.doubles( "O2_proliferation_saturation" );
    cell_defaults.parameters.o2_reference = cell_defaults.parameters.o2_proliferation_saturation;
    
    // set default motiltiy
    cell_defaults.phenotype.motility.is_motile = parameters.bools( "is_motile" );
    cell_defaults.functions.update_migration_bias = oxygen_taxis_motility;
    cell_defaults.phenotype.motility.persistence_time = parameters.doubles( "persistence_time" );
    cell_defaults.phenotype.motility.migration_bias = parameters.doubles( "migration_bias" );
    cell_defaults.phenotype.motility.migration_speed = parameters.doubles( "migration_speeed" ); // 0.5
    
    // set default uptake and secretion
    // oxygen
    cell_defaults.phenotype.secretion.secretion_rates[oxygen_i] = parameters.doubles( "O2_secretion_rate_constant" );
    cell_defaults.phenotype.secretion.uptake_rates[oxygen_i] = parameters.doubles( "O2_uptake_rate_constant" );
    cell_defaults.phenotype.secretion.saturation_densities[oxygen_i] = parameters.doubles( "O2_saturation_density" );
    
    // VEGF
    cell_defaults.phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
    cell_defaults.phenotype.secretion.uptake_rates[VEGF_i] = parameters.doubles("VEGF_uptake_rate_constant");
    cell_defaults.phenotype.secretion.saturation_densities[VEGF_i] = parameters.doubles("VEGF_saturation_density");
    
    // VEGF_Secretion and Vascular Death Function as Cell Rule
    // cell_defaults.functions.custom_cell_rule = VEGF_secretion_and_vascular_death_function;
    
    
    std::string tumorphenotype = parameters.strings("tumorphenotype");
    
    if (tumorphenotype == "model0")
    {
        cell_defaults.functions.update_phenotype = tumor_cell_phenotype0;
    }
    
    if (tumorphenotype == "model1")
    {
        cell_defaults.functions.update_phenotype = tumor_cell_phenotype1;
    }
    
    if (tumorphenotype == "model2")
    {
        cell_defaults.functions.update_phenotype = tumor_cell_phenotype2;
    }
    
    if (tumorphenotype == "model2a")
    {
        cell_defaults.functions.update_phenotype = tumor_cell_phenotype2a;
    }
    
    if (tumorphenotype == "model3")
    {
        cell_defaults.functions.update_phenotype = tumor_cell_phenotype3;
    }
    
    if (tumorphenotype == "model4")
    {
        cell_defaults.functions.update_phenotype = tumor_cell_phenotype4;
    }
    
    cell_defaults.name = "cancer cell";
    cell_defaults.type = 0;
    
    // add custom data
    
    std::vector<double> genes = { 1.0, 0.0 }; // RFP, GFP
    std::vector<double> proteins = {1.0, 0.0 }; // RFP, GFP;
    
    double default_degradation_rateRFP = parameters.doubles("default_degradation_rateRFP");
    double default_degradation_rateGFP = parameters.doubles("default_degradation_rateGFP");
    // 4.8e-4; // 24 hour half-life
    // 0.0077; // 90 minute half-life
    // 0.019; // 90% degrades in 120 minutes
    
    std::vector<double> degradation_rates = { default_degradation_rateRFP , default_degradation_rateGFP };
    
    double default_production_rateRFP = parameters.doubles("default_production_rateRFP");
    double default_production_rateGFP = parameters.doubles("default_production_rateGFP");    // 24 hour half ramp-up
    // 0.0019; // 6 hours to reach 50%
    // 0.0068; // 1.7 hours to reach 50%
    // 0.23; // 10 minutes to reach 90%
    
    std::vector<double> creation_rates = { default_production_rateRFP , default_production_rateGFP };
    
    cell_defaults.custom_data.add_vector_variable( "genes" , "dimensionless", genes );
    
    cell_defaults.custom_data.add_vector_variable( "proteins" , "dimensionless", proteins );
    cell_defaults.custom_data.add_vector_variable( "creation_rates" , "1/min" , creation_rates );
    cell_defaults.custom_data.add_vector_variable( "degradation_rates" , "1/min" , degradation_rates );
    
    cell_defaults.custom_data.add_variable( "hypoxic memory" , "min" , 0.0 );
    
    
    std::vector<double> color = {255, 255, 255};
    cell_defaults.custom_data.add_vector_variable( "nuclear_color" , "dimensionless", color );
    cell_defaults.custom_data.add_vector_variable( "cytoplasmic_color" , "dimensionless", color );
    
    return;
}

void setup_microenvironment( void )
{
    // set domain parameters
    
    default_microenvironment_options.X_range = {-1000, 1000};
    default_microenvironment_options.Y_range = {-1000, 1000};
    default_microenvironment_options.simulate_2D = true;
    
    // for 3a and 4
    default_microenvironment_options.X_range = {-1500, 1500};
    default_microenvironment_options.Y_range = {-1500, 1500};
    
    
    // gradients needed for this example
    
    default_microenvironment_options.calculate_gradients = true;
    
    // let BioFVM use oxygen as the default
    
    default_microenvironment_options.use_oxygen_as_first_field = true;
    
    // set Dirichlet conditions
    
    default_microenvironment_options.outer_Dirichlet_conditions = false;
    // default_microenvironment_options.Dirichlet_activation_vector[2] = false;
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("O2_Dirichlet_Condition");
    default_microenvironment_options.Dirichlet_activation_vector[0] = false;
    // default_microenvironment_options.Dirichlet_condition_vector[2] = 0.0;
    
    
    // add ECM
    
    Microenvironment* pME = get_default_microenvironment();
    
    /*     microenvironment.diffusion_coefficients[0]=10000000000000000;
     microenvironment.decay_rates[0]=10; */
    
    
    int ECM_i = pME->find_density_index( "ECM" );
    if( ECM_i < 0 )
    {
        std::cout << "Adding ECM to the microenvironment ... " << std::endl;
        pME->add_density( "ECM", "dimensionless" , 0.0 , 0.0 );
        ECM_i = pME->find_density_index( "ECM" );
        
        default_microenvironment_options.Dirichlet_condition_vector[ECM_i] = 1;
        default_microenvironment_options.Dirichlet_activation_vector[ECM_i] = false;
    }
    
    int VEGF_i = pME->find_density_index( "VEGF" );
    if( VEGF_i < 0 )
    {
        // 5.8 × 10−11 m2 s−1. // https://www.nature.com/articles/nprot.2012.051
        
        // decay 72h half-life : http://walter.deback.net/old/media/KohnLuque_etal_PhysBiol_2013.pdf
        // ECM binding: 1.5e-3 s-1 ~ 0.09 min-1
        
        std::cout << "Adding VEGF to the microenvironment ... " << std::endl;
        
        pME->add_density( "VEGF", "dimensionless" , 3.5e3 , 0.09 );
        VEGF_i = pME->find_density_index( "VEGF" );
        
        default_microenvironment_options.Dirichlet_condition_vector[VEGF_i] = 0.0;
        default_microenvironment_options.Dirichlet_activation_vector[VEGF_i] = false;
    }
    
    
    /*
     // add "repressor" signal
     
     pME->add_density( "Repressor" , 6.4e4 , .1 ); // 80 micron length scale
     // turn off Dirichlet for this factor
     default_microenvironment_options.Dirichlet_activation_vector[2] = false;
     */
    
    
    
    initialize_microenvironment();
    
    coarse_vasculature_setup(); // ANGIO
    
    
    // vascularization far initialization
    for(int i = 0; i<coarse_vasculature.vascular_densities.size();i++)
    {
        //   std::cout<<coarse_vasculature.vascular_densities[i].functional<<std::endl; 22500
        double tumor_radius = parameters.doubles("tumor_radius");
        double XPos=coarse_vasculature.mesh.voxels[i].center[0];
        double YPos=coarse_vasculature.mesh.voxels[i].center[1];
        double Center = 0;
        double Distance = XPos*XPos + YPos*YPos - Center;
        
        bool& A = pME->is_dirichlet_node(i);
        if (A == true)
        {
            pME->remove_dirichlet_node(i);
        }
        
        //std::cout<<A<<std::endl;
        if ( Distance < tumor_radius*tumor_radius)
        {
            coarse_vasculature.vascular_densities[i].functional = 0;
        }
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
    double tumor_confluence = parameters.doubles("tumor_confluence");
    
    double cell_spacing = tumor_confluence * 2.0 * cell_radius;
    
    double tumor_radius = parameters.doubles("tumor_radius");
    
    Cell* pCell = NULL;
    
    double x = 0.0;
    double x_outer = tumor_radius;
    double y = 0.0;
    
    double leader_fraction = parameters.doubles("leader_fraction"); // model 4
    // leader_fraction = 0.025; // model 4a  0.01 is a bit too small for this small starting group
    double follower_fraction = 1.0 - leader_fraction;
    
    /* pCell = create_cell(); // tumor cell
     pCell->assign_position( 0.0 , 0.0 , 0.0 ); */
    
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
            
            
            if( UniformRandom() <= follower_fraction)
            { pCell->type = 1; pCell->type_name = "Follower"; }
            
            pCell->custom_data[0] = NormalRandom( 1.0, 0.33 );
            if( pCell->custom_data[0] < 0.0 )
            { pCell->custom_data[0] = 0.0; }
            if( pCell->custom_data[0] > 2.0 )
            { pCell->custom_data[0] = .0; }
            
            
            if( fabs( y ) > 0.01 )
            {
                pCell = create_cell();  //tumor cell
                pCell->assign_position( x , -y , 0.0 );
                
                if( UniformRandom() <= follower_fraction)
                { pCell->type = 1; pCell->type_name = "Follower"; }
                
                pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
                if( pCell->custom_data[0] < 0.0 )
                { pCell->custom_data[0] = 0.0; }
                if( pCell->custom_data[0] > 2.0 )
                { pCell->custom_data[0] = .0; }
                
            }
            
            if( fabs( x ) > 0.01 )
            {
                pCell = create_cell(); // tumor cell
                pCell->assign_position( -x , y , 0.0 );
                
                if( UniformRandom() <= follower_fraction)
                { pCell->type = 1; pCell->type_name = "Follower"; }
                
                pCell->custom_data[0] = NormalRandom( 1.0, 0.25 );
                if( pCell->custom_data[0] < 0.0 )
                { pCell->custom_data[0] = 0.0; }
                if( pCell->custom_data[0] > 2.0 )
                { pCell->custom_data[0] = .0; }
                
                if( fabs( y ) > 0.01 )
                {
                    pCell = create_cell(); // tumor cell
                    pCell->assign_position( -x , -y , 0.0 );
                    
                    if( UniformRandom() <= follower_fraction)
                    { pCell->type = 1; pCell->type_name = "Follower"; }
                    
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
void tumor_cell_phenotype0( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 0
    // change color, but do nothing.
    
    static int genes_i = 0;
    static int proteins_i =1;
    static int creation_rates_i = 2;
    static int degradation_rates_i = 3;
    
    static int red_i = 0;
    static int green_i = 1;
    
    static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
    static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt); // call standard O2-based birth and death
    VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascularization
    
    // if cell is dead, don't bother with future phenotype changes.
    if( phenotype.death.dead == true )
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }
    
    // multiply proliferation rate by the oncoprotein
    
    //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    //static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    
    // set genes
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
    static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
    static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
    
    // permanent gene switch
    if( pO2 < FP_hypoxic_switch )
    {
        pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
        pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
    }
    
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
    }
    else
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
    return;
}


void tumor_cell_phenotype1( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 1
    // if pO2 < value, then motile, less proliferation (instantaneous)
    // if pO2 > value, then not motile, normal proliferation
    if( PhysiCell_globals.current_time < 720000 )
    {
        static int genes_i = 0;
        static int proteins_i =1;
        static int creation_rates_i = 2;
        static int degradation_rates_i = 3;
    
        static int red_i = 0;
        static int green_i = 1;
    
        static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
        static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
        static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
        update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
        VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
    
        // if cell is dead, don't bother with future phenotype changes.
        if( phenotype.death.dead == true )
        {
            pCell->functions.update_phenotype = NULL;
        return;
        }
    
        // multiply proliferation rate by the oncoprotein
        
        //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
        // int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
        
        //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
        
        // set genes
    
        double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
        static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
        static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
        
        // permanent gene switch
        if( pO2 < FP_hypoxic_switch )
        {
            pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
            pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
        }
        
        if( pO2 < phenotype_hypoxic_switch )
        {
            phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
        }
        else
        {
            phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
        // model 1
        // if pO2 < value, then motile, less proliferation (instantaneous)
        // if pO2 > value, then not motile, normal proliferation
        
        if( pO2 < phenotype_hypoxic_switch )
        {
            phenotype.motility.is_motile = true;
            // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
            phenotype.motility.migration_speed = parameters.doubles("migration_speed");
            // migration_bias = parameters.doubles("migration_bias");
        }
        else
        {
            phenotype.motility.is_motile = false;
        }
        
        return;
    }
    else
    {
        static int genes_i = 0;
        static int proteins_i =1;
        static int creation_rates_i = 2;
        static int degradation_rates_i = 3;
        
        static int red_i = 0;
        static int green_i = 1;
        
        static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
        
        static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
        static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
        
        update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
        VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
        
        int Ki67_negative_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_negative );
        int Ki67_positive_premitotic_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_positive_premitotic );
        int Ki67_positive_postmitotic_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_positive_postmitotic );
        
        phenotype.cycle.data.transition_rate( Ki67_negative_index ,Ki67_positive_premitotic_index ) = 0;  //0.004605;
        phenotype.cycle.data.transition_rate( Ki67_positive_premitotic_index ,Ki67_positive_postmitotic_index ) = 0;  //0.001282;
        phenotype.cycle.data.transition_rate( Ki67_positive_postmitotic_index, Ki67_negative_index ) = 0;  //0.0066667;
        
        int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
        
        phenotype.death.rates[apoptosis_index] = 0; //0.0000387
        
        // if cell is dead, don't bother with future phenotype changes.
        if( phenotype.death.dead == true )
        {
            pCell->functions.update_phenotype = NULL;
            return;
        }
        
        // multiply proliferation rate by the oncoprotein
        
        //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
        // int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
        
        //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
        
        // set genes
        
        double pO2 = (pCell->nearest_density_vector())[oxygen_i];
        
        static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
        static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
        
        // permanent gene switch
        if( pO2 < FP_hypoxic_switch )
        {
            pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
            pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
        }
        
        if( pO2 < phenotype_hypoxic_switch )
        {
            phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
        }
        else
        {
            phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
        
        // model 1
        // if pO2 < value, then motile, less proliferation (instantaneous)
        // if pO2 > value, then not motile, normal proliferation
        
        if( pO2 < phenotype_hypoxic_switch )
        {
            phenotype.motility.is_motile = false;
            // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
            phenotype.motility.migration_speed = parameters.doubles("migration_speed");
            // migration_bias = parameters.doubles("migration_bias");
        }
        else
        {
            phenotype.motility.is_motile = false;
        }
        
        return;
    }
}


void tumor_cell_phenotype2( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 2
    // if green gene on, then motile, less proliferation (permanent change)
    
    static int genes_i = 0;
    static int proteins_i =1;
    static int creation_rates_i = 2;
    static int degradation_rates_i = 3;
    
    static int red_i = 0;
    static int green_i = 1;
    
    static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
    static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
    VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
    
    // if cell is dead, don't bother with future phenotype changes.
    if( phenotype.death.dead == true )
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }
    
    // multiply proliferation rate by the oncoprotein
    
    //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    //static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    
    // set genes
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
    static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
    static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
    
    // permanent gene switch
    if( pO2 < FP_hypoxic_switch )
    {
        pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
        pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
    }
    
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
    }
    else
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
    // set phenotype in response to temporary or permanent changes
    
    // model 2
    // if green gene on, then motile, less proliferation (permanent change)
    
    if( pCell->custom_data.vector_variables[genes_i].value[green_i] > 0.1 )
    {
        phenotype.motility.is_motile = true;
        //        phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
        phenotype.motility.migration_speed = parameters.doubles("migration_speed");
        // migration_bias = parameters.doubles("migration_bias");
    }
    
    return;
}


void tumor_cell_phenotype2a( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 2a
    // if green gene on, then motile, proliferation unchanged (permanent change)
    // 2a: decrease adhesion
    
    static int genes_i = 0;
    static int proteins_i =1;
    static int creation_rates_i = 2;
    static int degradation_rates_i = 3;
    
    static int red_i = 0;
    static int green_i = 1;
    
    static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
    static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
    VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
    
    // if cell is dead, don't bother with future phenotype changes.
    if( phenotype.death.dead == true )
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }
    
    // multiply proliferation rate by the oncoprotein
    
    //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    //static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    
    // set genes
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
    static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
    static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
    
    // permanent gene switch
    if( pO2 < FP_hypoxic_switch )
    {
        pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
        pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
    }
    
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
    }
    else
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
    // set phenotype in response to temporary or permanent changes
    
    // model 2a
    // if green gene on, then motile, proliferation unchanged (permanent change)
    // 2a: decrease adhesion
    if( pCell->custom_data.vector_variables[genes_i].value[green_i] > 0.1 )
    {
        phenotype.motility.is_motile = true;
        // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
        phenotype.motility.migration_speed = parameters.doubles("migration_speed");
        // migration_bias = parameters.doubles("migration_bias");
        
        // phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("adhesion_distance");
    }
    return;
}


void tumor_cell_phenotype3( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 3
    // if pO2 < threshold, then set motile true, less proliferation
    // if pO2 > threshold and motile is true, probability rate*dt of turning motility false, restoring proliferation rate
    static int genes_i = 0;
    static int proteins_i =1;
    static int creation_rates_i = 2;
    static int degradation_rates_i = 3;
    
    static int red_i = 0;
    static int green_i = 1;
    
    static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
    static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
    VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
    
    // if cell is dead, don't bother with future phenotype changes.
    if( phenotype.death.dead == true )
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }
    
    // multiply proliferation rate by the oncoprotein
    
    //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    //static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    
    // set genes
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
    static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
    static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
    
    // permanent gene switch
    if( pO2 < FP_hypoxic_switch )
    {
        pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
        pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
    }
    
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
    }
    else
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
    // set phenotype in response to temporary or permanent changes
    
    
    // model 3
    // if pO2 < threshold, then set motile true, less proliferation
    // if pO2 > threshold and motile is true, probability rate*dt of turning motility false, restoring proliferation rate
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.motility.is_motile = true;
        phenotype.motility.migration_speed = parameters.doubles("migration_speed");
        // migration_bias = parameters.doubles("migration_bias");
        
        // phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("adhesion_distance");
    }
    else
    {
        static double persistence_time = parameters.doubles("persistence_time"); // 6 hours
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
    }
    return;
}


void tumor_cell_phenotype3a( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 3a
    // if pO2 < threshold, then set motile true, less proliferation
    // if pO2 > threshold and motile is true, probability rate*dt of turning motility false, restoring proliferation rate
    static int genes_i = 0;
    static int proteins_i =1;
    static int creation_rates_i = 2;
    static int degradation_rates_i = 3;
    
    static int red_i = 0;
    static int green_i = 1;
    
    static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
    static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
    VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
    
    // if cell is dead, don't bother with future phenotype changes.
    if( phenotype.death.dead == true )
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }
    
    // multiply proliferation rate by the oncoprotein
    
    //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    //static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    
    // set genes
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
    static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
    static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
    
    // permanent gene switch
    if( pO2 < FP_hypoxic_switch )
    {
        pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
        pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
    }
    
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
    }
    else
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
    // set phenotype in response to temporary or permanent changes
    
    
    // model 3a
    // if pO2 < threshold, then set motile true, less proliferation
    // if pO2 > threshold and motile is true, probability rate*dt of turning motility false, restoring proliferation rate
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.motility.is_motile = true;
        phenotype.motility.migration_speed = parameters.doubles("migration_speed");
        // migration_bias = parameters.doubles("migration_bias");
        
        // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
    }
    else
    {
        static double persistence_time = parameters.doubles("persistence_time"); //
        // 5760; // 4 days
        // 3600; // 60 hours
        // 360.0; // 6 hours // too short!
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
    }
    
    return;
}


void tumor_cell_phenotype4( Cell* pCell, Phenotype& phenotype, double dt )
{
    // model 4 and 4a
    // same as #3a, but only allowed in 10% of cells (permanent leader fraction)
    static int genes_i = 0;
    static int proteins_i =1;
    static int creation_rates_i = 2;
    static int degradation_rates_i = 3;
    
    static int red_i = 0;
    static int green_i = 1;
    
    static int hypoxic_memory_i = pCell->custom_data.find_variable_index( "hypoxic memory" );
    
    static int oxygen_i = get_default_microenvironment()->find_density_index( "oxygen" );
    static int VEGF_i = get_default_microenvironment()->find_density_index( "VEGF" );
    
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
    VEGF_secretion_and_vascular_death_function(pCell,phenotype,dt); // do vascular stuff
    
    // if cell is dead, don't bother with future phenotype changes.
    if( phenotype.death.dead == true )
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }
    
    // multiply proliferation rate by the oncoprotein
    
    //static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    //static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    //    phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    
    // set genes
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i];
    
    static double FP_hypoxic_switch = parameters.doubles("FP_hypoxic_switch");
    static double phenotype_hypoxic_switch = parameters.doubles("phenotype_hypoxic_switch");  //
    
    // permanent gene switch
    if( pO2 < FP_hypoxic_switch )
    {
        pCell->custom_data.vector_variables[genes_i].value[red_i] = 0.0;
        pCell->custom_data.vector_variables[genes_i].value[green_i] = 1.0;
    }
    
    if( pO2 < phenotype_hypoxic_switch )
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_hypoxia" );
    }
    else
    {
        phenotype.secretion.secretion_rates[VEGF_i] = parameters.doubles( "VEGF_secretion_rate_constant_normoxia" );
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
    
    // set phenotype in response to temporary or permanent changes
    
    
    // model 4 and 4a
    // same as #3a, but only allowed in 10% of cells (permanent leader fraction)
    if( pO2 < phenotype_hypoxic_switch && pCell->type != 1 )
    {
        phenotype.motility.is_motile = true;
        phenotype.motility.migration_speed = parameters.doubles("migration_speed");
        // migration_bias = parameters.doubles("migration_bias");
        // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
    }
    else
    {
        static double persistence_time = parameters.doubles("persistence_time"); //
        // 5760; // 4 days
        // 3600; // 60 hours
        // 360.0; // 6 hours // too short!
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
    
    // Necrotic - Purple
    if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling ||
       pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed ||
       pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
    {
        output[0] = "rgb(138,43,226)";
        output[2] = "rgb(148,0,211)";
        
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


void VEGF_secretion_and_vascular_death_function(Cell* pCell, Phenotype& phenotype, double dt )
{
    
    // VEGF SECRETION CODE
    /*
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
     
     // END VEGF SECRETION CODE */
    
    // VASCULAR DEATH CODE
    
    double vascular_degradation_rate_per_cell = parameters.doubles( "vascular_degradation_rate_per_cell" );
    
    coarse_vasculature( pCell ).functional = coarse_vasculature( pCell ).functional/(1+dt*vascular_degradation_rate_per_cell);
    
    // END VASCULAR DEATH CODE
    
    //    V - (V-1) = -dt*lambda*V
    //    V(1+dt*lambda) = (V-1)
    //    V=(V-1)/(1+dt*lambda)
    
    
}
