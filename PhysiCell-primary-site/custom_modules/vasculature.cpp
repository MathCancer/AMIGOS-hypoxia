//
//  vascularization.cpp
//  
//
//  Created by John Metzcar on 4/10/18.
//

#include "./vasculature.h"
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"




using namespace BioFVM;
using namespace PhysiCell;

#ifndef __vasculature_h__
#define __vasculature_h__

namespace PhysiCell{
    
Vascular_Options default_vascular_options;
Coarse_Vasculature coarse_vasculature;

Vascular_Options::Vascular_Options()
{
    vascular_mesh_multiplier = 1; //1;
    base_vascular_extension_rate = 0.0035/60.0;
    vascular_birth_rate = 1.0/18.0/6;
    max_vascular_density = 1.0; 
    vascular_death_rate = 1.0/18/60;  // 0
    vascular_proliferation_threshold = 0.001;  // 1
    vascular_proliferation_saturation= 0.5;   // 2
    vascular_chemotaxis_threshold = 0.001;   // 3
    vascular_chemotaxis_saturation = 0.5;   // 4
    effective_vascular_cutoff_threshold = 1.0E-8;  // 5
	
    
    angiogenesis_dt = 60.0;
    
    blood_substrate_densities.resize( 1 , 1.0 );
    tissue_far_field_substrate_densities.resize( 1 , 1.0 );
    
    blood_oxygen_tension = 38.0;
    tissue_far_field_oxygen_tension = 38.0;
    
//    degradation_rate_per_cell = 1e-6; // 1e-5; // How will we get cells to kill the vasculature? What can we use that is related to the cell presence/behavior and how to we access it?
    
    
    return;
}


void Vascular_Options::sync_to_BioFVM( void )
{
    // make sure it has the right number of densities.
    // set them equal to the boundary conditions for BioFVM
    
    blood_substrate_densities = default_microenvironment_options.Dirichlet_condition_vector;
    tissue_far_field_substrate_densities = default_microenvironment_options.Dirichlet_condition_vector;
    
    return;
}

Vascular_Densities::Vascular_Densities()
{
    functional = 0.5;
    total = 1.0;
    secretion_rate = 10.0; // could be too high
	target_O2 = 38.0;
	target_ECM = 0;
	target_VEGF = 0;
	target_vector = { target_O2,target_ECM,target_VEGF };
    vascular_extension_rate = 1.0;
//    vascular_birth_rate = 1.0/18.0;
    
    // How can I use the options more?
    
    return;
}

Coarse_Vasculature::Coarse_Vasculature()
{
    mesh.resize(1,1,1) ;
    vascular_densities.resize(1);
    
    blood_substrate_densities = default_vascular_options.blood_substrate_densities;
    
    
    VEGF.resize(1);
    net_vascular_density_fluxes.resize(1);
    
    
    pMicroenvironment = NULL;
    
    return;
}

void coarse_vasculature_setup( void )
{
    // sync the options to BioFVM
    default_vascular_options.sync_to_BioFVM();

    // use the custom bulk source/sink functions
    Microenvironment* pME = get_default_microenvironment();
    pME->bulk_supply_rate_function = vascular_supply_function;
    pME->bulk_supply_target_densities_function = vascular_target_function;
    pME->bulk_uptake_rate_function = vascular_uptake_function;

    // USER EDITS TO default_vascular_options GO HERE!!! (Why - becasue it is in the set up!) Duh!

    default_vascular_options.blood_oxygen_tension = 38;
    default_vascular_options.tissue_far_field_oxygen_tension = 38.0;
//
//    // END OF USER EDITS TO default_vascular_options
//
//    // sync the environment to BioFVM
//
    coarse_vasculature.sync_to_BioFVM();
//    coarse_vasculature.mesh.size();
//
//    // now, set bulk source and sink functions
//
    std::cout << "need to set up source functions!!!" << std::endl;
    
//    for(int i = 0; i<coarse_vasculature.vascular_densities.size();i++)
//    {
//        std::cout<<coarse_vasculature.vascular_densities[i].functional<<std::endl;
//    }
//
    return;
}

    
void Coarse_Vasculature::sync_to_BioFVM( void )
{
    // first, resize the mesh
    
    if( pMicroenvironment == NULL )
    { pMicroenvironment = get_default_microenvironment(); }
    
    int Xnodes = pMicroenvironment->mesh.x_coordinates.size();
    int Ynodes = pMicroenvironment->mesh.y_coordinates.size();
    int Znodes = pMicroenvironment->mesh.z_coordinates.size();
    
    if( Xnodes > 1 )
    { Xnodes /= default_vascular_options.vascular_mesh_multiplier; }
    if( Ynodes > 1 )
    { Ynodes /= default_vascular_options.vascular_mesh_multiplier; }
    if( Znodes > 1 )
    { Znodes /= default_vascular_options.vascular_mesh_multiplier; }
    
    // next, make sure the microenvironment has oxygen
    
/*     int oxygen_i = pMicroenvironment->find_density_index( "oxygen" );
    if( oxygen_i < 0 )
    {
        std::cout << "Adding oxygen to the microenvironment ... " << std::endl;
        
        if( default_microenvironment_options.use_oxygen_as_first_field == true )
        { pMicroenvironment->set_density( 0 , "oxygen", "mmHg" , 1e5 , 0.1 ); }
        else
        { pMicroenvironment->add_density( "oxygen", "mmHg" , 1e5 , 0.1 ); }
        oxygen_i = pMicroenvironment->find_density_index( "oxygen" );
        
//        default_microenvironment_options.Dirichlet_condition_vector[oxygen_i] = default_vascular_options.tissue_far_field_oxygen_tension;
//        default_microenvironment_options.Dirichlet_activation_vector[oxygen_i] = true;
    } */
    
    // next, make sure the microenvironment has VEGF
/*     
    int VEGF_i = pMicroenvironment->find_density_index( "VEGF" );
    if( VEGF_i < 0 )
    {
        // 5.8 × 10−11 m2 s−1. // https://www.nature.com/articles/nprot.2012.051
        
        // decay 72h half-life : http://walter.deback.net/old/media/KohnLuque_etal_PhysBiol_2013.pdf
        // ECM binding: 1.5e-3 s-1 ~ 0.09 min-1
        
        std::cout << "Adding VEGF to the microenvironment ... " << std::endl;
        
        pMicroenvironment->add_density( "VEGF", "dimensionless" , 3.5e3 , 0.09 );
        VEGF_i = pMicroenvironment->find_density_index( "VEGF" );
        
       default_microenvironment_options.Dirichlet_condition_vector[VEGF_i] = 0.0;
       default_microenvironment_options.Dirichlet_activation_vector[VEGF_i] = false;
    }
     */
    // next, resize the vascular mesh
    
    mesh.resize( default_microenvironment_options.X_range[0] , default_microenvironment_options.X_range[1] ,
                default_microenvironment_options.Y_range[0] , default_microenvironment_options.Y_range[1] ,
                default_microenvironment_options.Z_range[0] , default_microenvironment_options.Z_range[1] ,
                Xnodes, Ynodes, Znodes );
    mesh.units = default_microenvironment_options.spatial_units;
    
    // set the substrate densities to the correct values
    
    blood_substrate_densities = default_vascular_options.blood_substrate_densities;
    
    // next, make sure the vascular densities are of the right size
    
    vascular_densities.resize( mesh.voxels.size() );
    
    // now, make sure that the angiogenesis helper variables have the right size
    
    VEGF.resize( mesh.voxels.size() );
    net_vascular_density_fluxes.resize( mesh.voxels.size() );
    
    return;
}

Vascular_Densities& Coarse_Vasculature::operator()( std::vector<double> position )
{ return vascular_densities[ mesh.nearest_voxel_index( position ) ]; }

Vascular_Densities& Coarse_Vasculature::operator()( int n )
{ return vascular_densities[n]; }

Vascular_Densities& Coarse_Vasculature::operator()( int i, int j, int k )
{ return vascular_densities[ mesh.voxel_index(i,j,k) ]; }

Vascular_Densities& Coarse_Vasculature::operator()( Cell* pCell )
{ return this->operator()( pCell->position ); }

void Coarse_Vasculature::compute_coarse_VEGF( void )
{
    // if we're not synced to a microenvironment, then exit out
    if( pMicroenvironment == NULL )
    { return; }
    return;
}


void update_coarse_vasculature( double dt )
{
    
    // rho_v = .... see the CB code ...
    
//    if( pMicroenvironment == NULL )
//    { pMicroenvironment = get_default_microenvironment(); }
    
    Microenvironment* pMicroenvironment = get_default_microenvironment();
    update_vascular_extension_rate ( pMicroenvironment );
    update_vascular_population ( pMicroenvironment, dt );
    flux_vascular_density ( pMicroenvironment, dt )  ;
    
    //when will this get called???
    
    return;
}
    
void update_vascular_extension_rate ( Microenvironment* pMicroenvironment )  // Will move to Vasculature class
{
    
    extern Vascular_Options default_vascular_options;
    extern Coarse_Vasculature coarse_vasculature;
    
    static int VEGF_i = pMicroenvironment->find_density_index( "VEGF" );
    
    for ( int voxel_index = 0; voxel_index < coarse_vasculature.vascular_densities.size(); voxel_index++)
    {
         coarse_vasculature.vascular_densities[voxel_index].vascular_extension_rate = default_vascular_options.base_vascular_extension_rate*fmin(1 , fmax( 0, ( pMicroenvironment->density_vector(voxel_index)[VEGF_i]- default_vascular_options.vascular_chemotaxis_threshold)/(default_vascular_options.vascular_chemotaxis_saturation - default_vascular_options.vascular_chemotaxis_threshold)));
    }
    return;
}

void update_vascular_population ( Microenvironment* pMicroenvironment, double dt )
{
    
    double b, d;
    static int VEGF_i = pMicroenvironment->find_density_index( "VEGF" );
    extern Vascular_Options default_vascular_options;
    extern Coarse_Vasculature coarse_vasculature;
    
    
    
    for(int i = 0; i<coarse_vasculature.vascular_densities.size(); i++) // i ends up being the voxel index
    {
    
        b = default_vascular_options.vascular_birth_rate
            * fmax(0, ( 1 - coarse_vasculature.vascular_densities[i].functional / default_vascular_options.max_vascular_density ))
            * fmin(1 , fmax( 0, (pMicroenvironment->density_vector(i)[VEGF_i] - default_vascular_options.vascular_proliferation_threshold)/
            (default_vascular_options.vascular_proliferation_saturation - default_vascular_options.vascular_proliferation_threshold)));
    
        d = 0; //default_vascular_options.other_properties[0]*(PUT in SOMETHING RELATED TO CELLS HERE)/phenotypes_vector[0].max_cells;
    
        coarse_vasculature.vascular_densities[i].functional = coarse_vasculature.vascular_densities[i].functional / (1 - (b-d)*dt);
    
        if( coarse_vasculature.vascular_densities[i].functional < default_vascular_options.effective_vascular_cutoff_threshold)
        {
            coarse_vasculature.vascular_densities[i].functional = 0;
        }
    }
    return;
}

void flux_vascular_density ( Microenvironment* pMicroenvironment, double dt )  
{
    
    
//    std::vector<double> change_in_live_cell_population;
//
//    change_in_live_cell_population.assign(tissue_mesh.voxel_links.size(), 0.0);
    
    // Calculate fluxes - need it to look up, down, left, and right. How do I do that? Does BioFVM now its neighbors? What about just a straight vector instead of an array?  How about I make a links vector?
    
    // Just use the existing "inks" and calculate the flux into (or out of as the case might be each voxel at once by iterating thorugh conneted voxels list - since the rules require greater a to move, it hsould be one way always.... thie will require new code but shoudl run on older versions of BioFVM.
    
    extern Vascular_Options default_vascular_options;
    extern Coarse_Vasculature coarse_vasculature;
    
    int VEGF_i = pMicroenvironment->find_density_index( "VEGF" );
    
    double area = 20 * 20; // Exchange surface area of a 20 by 20 voxel - where is this actually stored????
    
    for( int voxel_index = 0; voxel_index < pMicroenvironment->mesh.connected_voxel_indices.size(); voxel_index++)
    {
        double a_i = pMicroenvironment->density_vector(voxel_index)[VEGF_i];
        double vascular_density_i = coarse_vasculature.vascular_densities[voxel_index].functional;
        double vascular_extention_rate_i = coarse_vasculature.vascular_densities[voxel_index].vascular_extension_rate;
//        std::cout<<coarse_vasculature.net_vascular_density_fluxes.size()<<std::endl;
        for( int nei_index = 0; nei_index < pMicroenvironment->mesh.connected_voxel_indices[voxel_index].size(); ++nei_index)
        {
//            std::cout<<nei_index<<std::endl;
            double a_j = pMicroenvironment->density_vector(nei_index)[VEGF_i];
        
                coarse_vasculature.net_vascular_density_fluxes[voxel_index] = -fmax(0,(a_j - a_i))
                *vascular_extention_rate_i * vascular_density_i + fmax(0,a_i-a_j)
                *coarse_vasculature.vascular_densities[nei_index].vascular_extension_rate *coarse_vasculature.vascular_densities[nei_index].functional;
        
        }
        
    }
    
        for( int voxel_index = 0; voxel_index < pMicroenvironment->mesh.connected_voxel_indices.size(); voxel_index++)
        {
//            Voxel* pV1 = tissue_mesh.voxel_links[k].pVoxel1;  // Is there a better way to store/access there rather than just access then twice in a row????
//            Voxel* pV2 = tissue_mesh.voxel_links[k].pVoxel2;
//            int i = pV1->mesh_index;
//            int j = pV2->mesh_index;
    
            coarse_vasculature.vascular_densities[voxel_index].functional = coarse_vasculature.vascular_densities[voxel_index].functional + coarse_vasculature.net_vascular_density_fluxes[voxel_index]*area*dt;

        }
//    std::cout<<pMicroenvironment->mesh.moore_connected_voxel_indices.size()<<std::endl;
//    return;
//
//    for( int k=0; k<tissue_mesh.voxel_links.size(); k++)
//    {
//        Voxel* pV1 = tissue_mesh.voxel_links[k].pVoxel1;
//        Voxel* pV2 = tissue_mesh.voxel_links[k].pVoxel2;
//        double area = tissue_mesh.voxel_links[k].surface_area;
//        double dx = tissue_mesh.voxel_links[k].center_to_center_distance;
//        int i = pV1->mesh_index;  // "testing" voxel
//        int j = pV2->mesh_index;  // "Target" voxel
//        double J_ij;
//
//        // k is the Link index!
//
//        // 02.23.18 eliminating this temporatliy. Vasculature SHOULD already be a denisty
//        //        double density1 = voxel_population_vectors[i].live_cell_counts[1]/voxel_population_vectors[i].phenotypes_vector[1].max_cells;
//        //        double density2 = voxel_population_vectors[j].live_cell_counts[1]/voxel_population_vectors[j].phenotypes_vector[1].max_cells;
//        // END elimination
//
//        // FROM 11/2
//
//        //        int h_ij = tissue_mesh.heaviside_fn(substrate_vectors[j].substrate_quantity[1] - substrate_vectors[i].substrate_quantity[1]);
//        //
//        //        int h_ji = tissue_mesh.heaviside_fn(substrate_vectors[i].substrate_quantity[1] - substrate_vectors[j].substrate_quantity[1]);
//        //        Variant 1 - "diffusive rho_v" 11-2-17
//
//        //        J_ij = -fmax(0,(density1 - density2))*voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])*h_ij
//        //        + fmax(0,density2-density1)*voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1])*h_ji;
//
//        //        Variant 2 - non-difference rho_v, perhaps more advective 11-29-17
//
//        //        J_ij = -fmax(0,(density1))*voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])*h_ij
//        //        + fmax(0,density2)*voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1])*h_ji;
//
//        //        std::cout<<voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])<<std::endl;
//        // END 11/2
//
//        // FROM 11/7 - Why would this not work on a sphere?  - well if the diffusion equations aren't right ... like if you accidently get a reversed gradient ... duh.  Fix that first.
//
//        // Updated to include correct density calculation 1/6/18
//
//
//        // 02.23.18 eliminating this temporatliy. Vasculature SHOULD already be a denisty AND AF should also arelady be a density
//        //        J_ij = -fmax(0,(substrate_vectors[j].substrate_quantity[1]/tissue_mesh.voxels[j].volume - substrate_vectors[i].substrate_quantity[1]/tissue_mesh.voxels[i].volume))*voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])
//        //        *density1 + fmax(0,substrate_vectors[i].substrate_quantity[1]/tissue_mesh.voxels[i].volume-substrate_vectors[j].substrate_quantity[1]/tissue_mesh.voxels[j].volume)
//        //        *voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1])*density2;
//        // 02.23.18 END old old
//
//
//        // 02.23.18 updated flux code - reflecting the TRUTH that AF and vasculature are densities
//
//        J_ij = -fmax(0,(substrate_vectors[j].substrate_quantity[1] - substrate_vectors[i].substrate_quantity[1]))
//        *voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])
//        *voxel_population_vectors[i].live_cell_counts[1] + fmax(0,substrate_vectors[i].substrate_quantity[1]-substrate_vectors[j].substrate_quantity[1])
//        *voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1]) *voxel_population_vectors[j].live_cell_counts[1];
//
//
//        // END update
//
//        // END 11/7
//
//        // From 11/30/17  - will test later
//
//        // Variant 1 - "diffusive rho_v" (like from 11/7) but with angioF scaled by the value in the voxel versus the range of the substart AND taking a density difference instead of being just the denisty
//
//        // Varient 2 - "advective rho_v" like 11/7 but with angioF scaled by the value in the voxel versus the range of the substart
//
//        //        J_apoptotic_cells[i] = (density2 - density1)/dx;
//        //        J_necrotic_cells[i] = (density2 - density1)/dx;
//
//        change_in_live_cell_population[k] = -dt*J_ij*area;
//
//    }
//
//    for( int k=0; k<tissue_mesh.voxel_links.size(); k++)
//    {
//        Voxel* pV1 = tissue_mesh.voxel_links[k].pVoxel1;  // Is there a better way to store/access there rather than just access then twice in a row????
//        Voxel* pV2 = tissue_mesh.voxel_links[k].pVoxel2;
//        int i = pV1->mesh_index;
//        int j = pV2->mesh_index;
//
//        voxel_population_vectors[i].live_cell_counts[1] = voxel_population_vectors[i].live_cell_counts[1] - change_in_live_cell_population[k];
//        voxel_population_vectors[j].live_cell_counts[1] = voxel_population_vectors[j].live_cell_counts[1] + change_in_live_cell_population[k];
//    }
//
    
    return;
}
    
void vascular_supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{

        // use this syntax to access the jth substrate write_here
        // (*write_here)[j]
        // use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
        // microenvironment->density_vector(voxel_index)[j]

        extern Vascular_Options default_vascular_options;
        extern Coarse_Vasculature coarse_vasculature;

//        static std::vector<double> delivery_rate_vector( microenvironment->number_of_densities() , 1000.0 );
//        static bool setup_done = false;
//        if( setup_done == false )
//        {
//            for( int i=0; i < microenvironment->number_of_densities() ; i++ )
//            {
//                if( default_microenvironment_options.Dirichlet_activation_vector[i] == false )
//                { delivery_rate_vector[i] = 0.0; }
//            }
//            setup_done = true;
//        }
//
        // functional_vascular_density * source_rates * on_off



        for( int i=0 ; i < 1; i++ )
        {
            (*write_here)[i] = coarse_vasculature( microenvironment->mesh.voxels[voxel_index].center ).functional; // what is it doing?
            (*write_here)[i] *= coarse_vasculature( microenvironment->mesh.voxels[voxel_index].center ).secretion_rate ;
//            std::cout<<i<<std::endl;
        }
        return;


    return;
}
    
//void vascular_supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
//{
//    // use this syntax to access the jth substrate write_here
//    // (*write_here)[j]
//    // use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
//    // microenvironment->density_vector(voxel_index)[j]
//
//    extern Vascular_Options default_vascular_options;
//    extern Coarse_Vasculature coarse_vasculature;
//
//    /*
//     // DEBUG
//     for( int i=0 ; i < write_here->size() ; i++ )
//     { (*write_here)[i] = 0.0; }
//     return;
//     */
//
//    // figure out which substrate get delivered to/from the vasculature
//    static std::vector<double> delivery_rate_vector( microenvironment->number_of_densities() , 1000.0 );
//    static bool setup_done = false;
//    if( setup_done == false )
//    {
//        for( int i=0; i < microenvironment->number_of_densities() ; i++ )
//        {
//            if( default_microenvironment_options.Dirichlet_activation_vector[i] == false )
//            { delivery_rate_vector[i] = 0.0; }
//        }
//        setup_done = true;
//    }
//
//    // functional_vascular_density * source_rates * on_off
//
//
//
//    for( int i=0 ; i < write_here->size() ; i++ )
//    {
//        (*write_here)[i] = coarse_vasculature( microenvironment->mesh.voxels[voxel_index].center ).functional;
//        (*write_here)[i] *= delivery_rate_vector[i];
//    }
//    return;
//
//
//    //    *(write_here) = delivery_rate_vector; // S_i = delivery_rate
//    //    *(write_here) *= coarse_vasculature( microenvironment->mesh.voxels[voxel_index].center ).functional; // S_i = functional(i)*delivery_rate;
//
//    return;
//}

void vascular_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{

    //    // use this syntax to access the jth substrate write_here
    //    // (*write_here)[j]
    //    // use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
    //    // microenvironment->density_vector(voxel_index)[j]
    //
    extern Coarse_Vasculature coarse_vasculature;
    extern Vascular_Options default_vascular_options;
/* 
	
	double target_O2 = 38.0;
	double target_ECM = 0;
	double target_VEGF = 0;
	std::vector<double> target_vector = { target_O2,target_ECM,target_VEGF }; */


    for( int i=0 ; i < write_here->size() ; i++ )
    { (*write_here)[i] = coarse_vasculature( microenvironment->mesh.voxels[voxel_index].center ).target_vector[i] ; }
    return;
}
    
//void vascular_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
//{
//    // use this syntax to access the jth substrate write_here
//    // (*write_here)[j]
//    // use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
//    // microenvironment->density_vector(voxel_index)[j]
//
//    extern Coarse_Vasculature coarse_vasculature;
//    extern Vascular_Options default_vascular_options;
//
//    for( int i=0 ; i < write_here->size() ; i++ )
//    { (*write_here)[i] = coarse_vasculature.blood_substrate_densities[i]; }
//    return;
//    /*
//     for( int i=0 ; i < write_here->size() ; i++ )
//     { (*write_here)[i] = coarse_vasculature.blood_substrate_densities[i]; }
//     return;
//     */
//
//    //    (*write_here) = default_microenvironment_options.Dirichlet_condition_vector;
//
//    return;
//}


void vascular_uptake_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
    //    static std::vector<double> uptake_rate_vector( 1.0 , microenvironment->number_of_densities() );
    
    // DEBUG
//    for( int i=0 ; i < write_here->size() ; i++ )
//    { (*write_here)[i] = 0.0; }
//    return;
//
//
//    for( int i=0 ; i < write_here->size() ; i++ )
//    { (*write_here)[i] = 0.0; }
//
//    //    (*write_here) = uptake_rate_vector;
    return;
}
    
void write_vasculature_data_matlab( std::string filename )

{
    
    int number_of_data_entries = microenvironment.number_of_voxels();
    
    int size_of_each_datum = 6;
    
    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "Vascular_Data" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.
    
    for( int i=0; i < number_of_data_entries ; i++ )
        
    {
//        double temp1 = microenvironment.mesh.voxels[i].center[0];

        fwrite( (char*) &( coarse_vasculature.mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( coarse_vasculature.mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( coarse_vasculature.mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( coarse_vasculature.mesh.voxels[i].volume ) , sizeof(double) , 1 , fp );

        fwrite( (char*) &( coarse_vasculature.vascular_densities[i].functional), sizeof(double) , 1 , fp );
        
        fwrite( (char*) &( coarse_vasculature.vascular_densities[i].total), sizeof(double) , 1 , fp );
        
        // current voxel index of cell
        
    }
    
    
    
    fclose( fp );
    
    
    
    return;
    
}
    
    
}

#endif
