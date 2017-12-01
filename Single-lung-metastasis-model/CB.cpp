//
//  CB.cpp
//  
//
//  Created by John Metzcar on 8/30/17.
//
//

#include "./BioFVM.h"  // SHould I just include all the headers?  Why or why not?

//#include "BioFVM_vector.h"
//#include "BioFVM_mesh.h"
//#include "CB.h"

namespace BioFVM{

Phenotype::Phenotype()
 {
     name = "default"; // From previous Matlab prototyping.  Roughly based on PM's previous work/PhysiCell.
     
     // Do we need name here?  Does it really get stored as the actual name of the instantiation?  Like "Phenotype MCF-7;"???
     
     phenotype_ID = 0;
     
     time_units = "min";
     spatial_units ="micron";
     
     birth_rate = 1.0/(18.0*60);
     apoptotic_death_rate = 0.01*birth_rate;
     apoptotic_clearance_rate = 1.0/(60*8.6);
     necrotic_death_rate = 1.0/(60*1800.0);//0.01*birth_rate;//0.5; //%1/(1*24);%1/(3*24.0); % 1.0/(3 * 24.0 ); % 3 day survival ;
     necrotic_clearance_rate = 1.0/(60*60.0*24.0);
     motility = 135.0/60; // um^2/hr?
     spatial_mechanical_factor = 1.05; // Gives the percentage above maximum cell packing that cells can grow to (on the voxel level)
     spatial_proliferation_factor = 0.95; // Gives the percentage of maxium cell packing at which the cells start to spill into neighboring voxels;
     base_secretion_rate = 1;  // Base/max AF secretion rate (dimensionless)
     hypoxic_o2_threshold = 5.0/38.0; // Matches PC hypoxic threshold for breast cancer, oxygen well vasculartized equals 38 mmHg
     critical_o2_threshold = 2.5/38.0; // Matches PC critical for breast cancer, oxygen well vasculartized equals 38 mmHg
     cell_volume = (4/3)*3.1415926*10E3;  //
     max_cells = 216.0;  // 60*60*60 (voxel volume)/10*10*10 (cell volume)  all in microns
     other_properties.resize(0);
     other_properties_names.resize(0);
     
     std::vector<double> temp;
     temp.assign(3, 0.0);
     
     for(int i = 0; i<2; i++)
     {
         secretion_rate.push_back(temp);
         substrate_target.push_back(temp);
         uptake_rate.push_back(temp);
     }

     
     /* Maybe move these to a vector such that they can be easily added?  Or maybe add an unspecificed vector that can hodl more as required?*/
     
     return;
     
 }


Voxel_Population_Vector::Voxel_Population_Vector()
 {
     tissue_voxel_index = 0;
     pMicroenvironment = NULL;
     microenvironment_voxel_indices.resize(0);
     names_vector.resize(0);  // Will be moved out eventually
     phenotypes_vector.resize(0); // vector of each phenotype present in voxel
     live_cell_counts.resize(0);  // Vector of live cell counts for each population (on per voxel level)
     apoptotic_cell_counts.resize(0);  // Vector of apoptotic cell counts for each population (on per voxel basis) - Could see this going in as just a population instead of subpopulation as needed
     necrotic_cell_counts.resize(0);  // Vector of necrotic_cell_counts for each population (on per voxel basis) - Could see this going in as just a population instead of subpopulation as needed
     
     return;
 }
    
void Voxel_Population_Vector::set_population(int population_number, double live_cell_count, std::string name) 
 {
     names_vector[population_number] = name;
     live_cell_counts[population_number] = live_cell_count;
     
     return;
 }

void Voxel_Population_Vector::set_population(int population_number, double live_cell_count, double apoptotic_cell_count, double necrotic_cell_count, std::string name)
 {
     names_vector[population_number] = name;
     live_cell_counts[population_number] = live_cell_count;
     apoptotic_cell_counts[population_number] = apoptotic_cell_count;
     necrotic_cell_counts[population_number] = necrotic_cell_count;
     
     return;
 }
    
void Voxel_Population_Vector::add_population(std::string name, Phenotype phenotype)
 {
     names_vector.push_back(name);
     phenotypes_vector.push_back(phenotype);
     phenotype.phenotype_ID = phenotypes_vector.size()-1;
     live_cell_counts.push_back(0);
     apoptotic_cell_counts.push_back(0);
     necrotic_cell_counts.push_back(0);
     
     return;
 }

void Voxel_Population_Vector::pair_to_microenvironment( Microenvironment* pME )
 {
     pMicroenvironment = pME;
     
     return;
 }

std::vector<double> Voxel_Population_Vector::sample_microenvironment( void )  // Not currently using.  Not using the microenvironemnt class from BioFVM.  Same technique could be useful in bioTissueBox, but not implemented yet.
{
    std::vector<double> output;
    
    if( pMicroenvironment == NULL )
    {
        output.resize(0);
        return output;
    }
    
    output.resize( pMicroenvironment->number_of_densities() , 0.0 );
    
    double total_volume = 0.0;
    
    for( int i=0 ; i < microenvironment_voxel_indices.size() ; i++ )
    {
        int n = microenvironment_voxel_indices[i];
        // output += (*pMicroenvironment)(n); // get the substrate density vector at voxel j
        
        // output = output + density(n)*Volume(n);
        
        total_volume += pMicroenvironment->mesh.voxels[n].volume;
        // output = output + pMicroenvironment->mesh.voxels[n].volume * (*pMicroenvironment)(n);
        // // y = y + a*x
        // void axpy( std::vector<double>* y, double& a , std::vector<double>& x );
        axpy( &output , pMicroenvironment->mesh.voxels[n].volume , (*pMicroenvironment)(n) );  // Why does thie function and this class have access to this function?
        // Also and more importantly, what is a functor and what is going on in lines 603 and 604 of uE.cpp?  How do we know which density substrate vector it used?   
    }
    output /= ( 1e-16 + total_volume );
    return output;
    
}
    
double Voxel_Population_Vector::update_AF_secretion_rate ( double oxygen )
{
    double modified_rate;

    if( oxygen > phenotypes_vector[0].hypoxic_o2_threshold )
    {

        modified_rate = 0.0;}
    
    else if ( oxygen < phenotypes_vector[0].critical_o2_threshold)
    {

        modified_rate = phenotypes_vector[0].base_secretion_rate;}
    
    else
    {

        modified_rate = phenotypes_vector[0].base_secretion_rate*((phenotypes_vector[0].hypoxic_o2_threshold - oxygen)/(phenotypes_vector[0].hypoxic_o2_threshold - phenotypes_vector[0].critical_o2_threshold));
    }
    
    
    return modified_rate;
}
    
double Voxel_Population_Vector::update_vascular_creep_rate ( double AF )  // Will move to Vasculature class
{
    
    double modified_rate = phenotypes_vector[1].motility*fmin(1 , fmax( 0, (AF - phenotypes_vector[1].other_properties[3])/(phenotypes_vector[1].other_properties[4] - phenotypes_vector[1].other_properties[3])));
    
    return modified_rate;
}
    
void Voxel_Population_Vector::update_vascular_population ( double dt, double AF )
{
    
    double b, d;
    
    b = phenotypes_vector[1].birth_rate * ( 1 - live_cell_counts[1] / phenotypes_vector[1].max_cells )
        *fmin(1 , fmax( 0, (AF - phenotypes_vector[1].other_properties[1])/(phenotypes_vector[1].other_properties[2] - phenotypes_vector[1].other_properties[1])));
    
    d = phenotypes_vector[1].other_properties[0]*(live_cell_counts[0]+apoptotic_cell_counts[0]+necrotic_cell_counts[0])/phenotypes_vector[0].max_cells;
    
    live_cell_counts[1] = live_cell_counts[1] / (1 - (b-d)*dt);
    
    if( live_cell_counts[1] < phenotypes_vector[1].other_properties[5] )
    {
        live_cell_counts[1] = 0;
    }
    
    return;
}
    
void Voxel_Population_Vector::update_populations( double dt, double oxygen )  // 1) Pull out the appropriate values from the phenotype instants.
                                      // 2) Also need to be able to support more than one type of population.
 {
     int i = 0;  // NOTE WILL NEED MODIFIED TO DO MORE THAN 1 POPULATION!!
     double b, d, d_necrotic;
     
     b = phenotypes_vector[i].birth_rate * fmax(0 , (oxygen -  phenotypes_vector[i].hypoxic_o2_threshold )/(1-phenotypes_vector[i].hypoxic_o2_threshold)) * fmax(0,(1.0 -  live_cell_counts[i] /  ( phenotypes_vector[i].spatial_proliferation_factor * phenotypes_vector[i].max_cells)) );
          
     d_necrotic = phenotypes_vector[i].necrotic_death_rate * fmin( 1 , ((phenotypes_vector[i].hypoxic_o2_threshold - oxygen)/(phenotypes_vector[i].hypoxic_o2_threshold - phenotypes_vector[i].critical_o2_threshold)));
          
     if( d_necrotic < 0 )
     {
          d_necrotic = 0;
     }

     d = phenotypes_vector[i].apoptotic_death_rate + d_necrotic;
          
     live_cell_counts[i] = live_cell_counts[i] / (1 - (b-d)*dt);
     apoptotic_cell_counts[i] = (apoptotic_cell_counts[i]+dt*phenotypes_vector[i].apoptotic_death_rate*live_cell_counts[i]) / ( 1+phenotypes_vector[i].apoptotic_clearance_rate*dt );
     necrotic_cell_counts[i] = (necrotic_cell_counts[i]+dt*d_necrotic*live_cell_counts[i]) / ( 1+phenotypes_vector[i].necrotic_clearance_rate*dt );
     
     return;
 }

Substrate_Vector::Substrate_Vector()
{
    
    int tissue_mesh_index = 0;
    names_vector.resize(0);  // Will be moved out eventually
    units_vector.resize(0);  // Will be moved out eventually
    substrate_quantity.resize(0);
    substrate_decay_constant.resize(0);
    
    return;
}

void Substrate_Vector::add_substrate(std::string name, std::string unit)
{
    names_vector.push_back(name);
    units_vector.push_back(unit);
    substrate_quantity.push_back(0);
    substrate_decay_constant.push_back(0);
    
    return;
}

void Substrate_Vector::set_substrate(int substrate_number, double quantity, std::string name, std::string unit)
{
    
    names_vector[substrate_number] = name;
    units_vector[substrate_number] = unit;
    substrate_quantity[substrate_number] = quantity;
    
    return;
}

void Substrate_Vector::set_decay_constant (int substrate_number, double decay_constant)
{

    substrate_decay_constant[substrate_number] = decay_constant;
    
    return;
}
Substrate_Properties::Substrate_Properties()
{
    diffusion_coefficients.resize(0);
    
    return;
}
    
void Substrate_Properties::add_diffusion_coefficient(double diffusion_coefficient)
{
    diffusion_coefficients.push_back(diffusion_coefficient);
}
    
Tissue::Tissue()
 {

     General_Mesh tissue_mesh;
     Microenvironment chemical_microenvironment;
     //tissue_mesh.voxels.clear();  // IN the default constructor, the general mesh has one default voxel.  This effectively overrides that for the Tissue Mesh.  Not sure if that one voxel is required for somethign else.  WATCH FOR BUGs!!
     voxel_population_vectors.resize(0);
     vascular_densities.resize(0);
     max_vascular_densities.resize(0);
     

     
     return;
 }
    
void Tissue::create_population_vectors (void) // but then all the phenotypes need assigned.  Can do it this way but it might not be the best.
{
    for(int i = 0; i<tissue_mesh.voxels.size(); i++)
    {
        Voxel_Population_Vector population_vector;
        population_vector.tissue_voxel_index = tissue_mesh.voxels[i].mesh_index;
        population_vector.pMicroenvironment = &chemical_microenvironment;
        voxel_population_vectors.push_back( population_vector);
    }
    return;
}
    
void Tissue::create_substrate_vectors (void)
{
    for(int i = 0; i<tissue_mesh.voxels.size(); i++)
    {
        Substrate_Vector substrate_vector;
        substrate_vector.tissue_voxel_index = tissue_mesh.voxels[i].mesh_index;
        substrate_vectors.push_back( substrate_vector );
    }
    return;
}
    
void Tissue::create_substrate_properties_vectors (void)
{
    for(int i=0; i<tissue_mesh.voxel_links.size(); i++)
    {
        Substrate_Properties substrate_properties;
        substrate_properties_vectors.push_back(substrate_properties);
    }

    return;
}

void Tissue::make_tumor_spheroid_2D( double tumor_radius , double voxel_radius )  // Hexognal lattice  Could flag this is desired to make it work for either 2 or 3 D?  Currently set only for 2-D
    {
        
        double x_position = -tumor_radius;
        double y_position = -tumor_radius;
        double z_position = 0;//-tumor_radius;  // This to non-zero to have a 3 D mesh.
        double V = 4/3*3.141592654*pow(voxel_radius,3);
        
        double x_offset = 0;
        int row_count = 0;
        
        std::cout << "I'm starting" << std::endl;
        
        while( y_position <= tumor_radius + 1e-15 )
        {
            
            while( x_position <= tumor_radius + 1e-15 )
            {
                // drop a voxel at (x_position,y_position)
                
                // keep the voxel if it's inside the circle
                if( x_position * x_position + y_position*y_position + z_position*z_position<= tumor_radius*tumor_radius )
                {
                    //                std::cout << x_position << "," << y_position << "\t";
                    
                    tissue_mesh.add_voxel(x_position,y_position,z_position,V);
                    
                }
                
                // march along in x
                x_position += 2.0 * voxel_radius;
            }
            
            y_position += sqrt(3)*voxel_radius;
            row_count++;
            
            if( row_count % 2 == 1 )
            { x_offset = voxel_radius; }
            else
            { x_offset = 0.0; }
            
            x_position = -tumor_radius + x_offset;
        }
        
        std::cout << "Spheroid set up on hex lattice" << std::endl;
        
        return;
        
    }

void Tissue::link_tumor_spheroid( double test_distance, double shared_surface_area )  // Notice all the declarations that are immediately defined.  I wonder if that is the best practice.
{

    std::vector<std::vector<double>> used_indices;
    
    for(int i=0; i<tissue_mesh.voxels.size(); i++)
    {
        double x_coordinate_i = tissue_mesh.voxels[i].center[0];
        double y_coordinate_i =tissue_mesh.voxels[i].center[1];
        double z_coordinate_i =tissue_mesh.voxels[i].center[2];
        
        for(int j=i; j<tissue_mesh.voxels.size(); j++)
        {
            double x_coordinate_j = tissue_mesh.voxels[j].center[0];
            double y_coordinate_j = tissue_mesh.voxels[j].center[1];
            double z_coordinate_j = tissue_mesh.voxels[j].center[2];
            
            double distance_sqr = pow(x_coordinate_i - x_coordinate_j,2) + pow(y_coordinate_i-y_coordinate_j,2) + pow(z_coordinate_i-z_coordinate_j,2);
            double test_distance_sqr = pow(test_distance,2);
            
            if((distance_sqr<test_distance_sqr+1E-2)&&(i!=j))
            {
                double center_to_center_distance = sqrt(distance_sqr);
                std::vector<double> normal_vector_i_to_j;
                normal_vector_i_to_j.push_back((x_coordinate_j-x_coordinate_i)/(center_to_center_distance+1E-14));
                normal_vector_i_to_j.push_back((y_coordinate_j - y_coordinate_i)/(center_to_center_distance+1E-14));
                normal_vector_i_to_j.push_back((z_coordinate_j - z_coordinate_i)/(center_to_center_distance+1E-14));
//                std::cout<<(x_coordinate_j-x_coordinate_i + y_coordinate_j - y_coordinate_i)/center_to_center_distance<<std::endl;
                tissue_mesh.link_voxels(tissue_mesh.voxels[i].mesh_index, tissue_mesh.voxels[j].mesh_index, shared_surface_area, center_to_center_distance, normal_vector_i_to_j);
            }
            
        }
        
    }
    
    
    return;
}

void Tissue::update_substrates( double dt )  // Needs targets, secretions, and decays for all substrates, cell types, and cell fractions.  Those values should go into the phenotype vector for the voxel.
{
    for(int k=0; k<substrate_vectors.size(); k++)  // selects all individually stored substrate vectors (accesses all voxels)
    {

        for( int j=0; j<voxel_population_vectors[k].phenotypes_vector.size(); j++)  // Selects all phenotypes/cell types in a voxel
        {

            for(int i=0; i<substrate_vectors[k].substrate_quantity.size(); i++)  // Selects all the substrates and adds across phenotype subtypes (live, apop, nec)
            {

                //  Add sums for the subtypes (apoptotic, necrotic, etc).
                double temp_secretion = (voxel_population_vectors[k].phenotypes_vector[j].secretion_rate[i][0] * voxel_population_vectors[k].live_cell_counts[j] +
                                         voxel_population_vectors[k].phenotypes_vector[j].secretion_rate[i][1] * voxel_population_vectors[k].apoptotic_cell_counts[j] +
                                         voxel_population_vectors[k].phenotypes_vector[j].secretion_rate[i][2] * voxel_population_vectors[k].necrotic_cell_counts[j])/
                                         voxel_population_vectors[k].phenotypes_vector[j].max_cells;

                double temp_target = (voxel_population_vectors[k].phenotypes_vector[j].substrate_target[i][0] * voxel_population_vectors[k].live_cell_counts[j] +
                                      voxel_population_vectors[k].phenotypes_vector[j].substrate_target[i][1] * voxel_population_vectors[k].apoptotic_cell_counts[j] +
                                      voxel_population_vectors[k].phenotypes_vector[j].substrate_target[i][2] * voxel_population_vectors[k].necrotic_cell_counts[j])/
                                      voxel_population_vectors[k].phenotypes_vector[j].max_cells;

                double temp_uptake = (voxel_population_vectors[k].phenotypes_vector[j].uptake_rate[i][0] * voxel_population_vectors[k].live_cell_counts[j] +
                                      voxel_population_vectors[k].phenotypes_vector[j].uptake_rate[i][1] * voxel_population_vectors[k].apoptotic_cell_counts[j] +
                                      voxel_population_vectors[k].phenotypes_vector[j].uptake_rate[i][2] * voxel_population_vectors[k].necrotic_cell_counts[j])/
                                      voxel_population_vectors[k].phenotypes_vector[j].max_cells;

                substrate_vectors[k].substrate_quantity[i] = (substrate_vectors[k].substrate_quantity[i] + dt * temp_secretion * temp_target)/(1+dt*(temp_secretion+temp_uptake+substrate_vectors[k].substrate_decay_constant[i]));  // DOUBLE CHECK EQUATION!!!
            }
        }
    }
    return;
}
    
void Tissue::run_diffusion ( double dt, int substrate_index)  // Is there any reason to not just loop over the size of substrate vector?
{
    std::vector<double> change_in_substrate;
    change_in_substrate.assign(tissue_mesh.voxel_links.size(), 0.0);

    // Calculate fluxes
    
    for( int i=0; i<tissue_mesh.voxel_links.size(); i++)
    {
        Voxel* pV1 = tissue_mesh.voxel_links[i].pVoxel1;
        Voxel* pV2 = tissue_mesh.voxel_links[i].pVoxel2;
        double area = tissue_mesh.voxel_links[i].surface_area;
        double dx = tissue_mesh.voxel_links[i].center_to_center_distance;
        int voxel1_index = pV1->mesh_index;
        int voxel2_index = pV2->mesh_index;
        
        double density1 = (substrate_vectors[voxel1_index].substrate_quantity[substrate_index])/tissue_mesh.voxels[voxel1_index].volume;  
        double density2 = (substrate_vectors[voxel2_index].substrate_quantity[substrate_index])/tissue_mesh.voxels[voxel2_index].volume;
        
        change_in_substrate[i] = -dt*(density2 - density1)/dx*area*substrate_properties_vectors[i].diffusion_coefficients[0];
        
    }
    
    for( int i=0; i<tissue_mesh.voxel_links.size(); i++)
    {
        Voxel* pV1 = tissue_mesh.voxel_links[i].pVoxel1;  // Is there a better way to store/access there rather than just access then twice in a row????
        Voxel* pV2 = tissue_mesh.voxel_links[i].pVoxel2;
        int voxel1_index = pV1->mesh_index;
        int voxel2_index = pV2->mesh_index;
        
        substrate_vectors[voxel1_index].substrate_quantity[substrate_index] = substrate_vectors[voxel1_index].substrate_quantity[substrate_index] - change_in_substrate[i];
        substrate_vectors[voxel2_index].substrate_quantity[substrate_index] = substrate_vectors[voxel2_index].substrate_quantity[substrate_index] + change_in_substrate[i];
    }
    return;
}

void Tissue::flux_vascular_density ( double dt )  // Will be moved to vasculature class
{


    std::vector<double> change_in_live_cell_population;

    change_in_live_cell_population.assign(tissue_mesh.voxel_links.size(), 0.0);

    // Calculate fluxes

    for( int k=0; k<tissue_mesh.voxel_links.size(); k++)
    {
        Voxel* pV1 = tissue_mesh.voxel_links[k].pVoxel1;
        Voxel* pV2 = tissue_mesh.voxel_links[k].pVoxel2;
        double area = tissue_mesh.voxel_links[k].surface_area;
        double dx = tissue_mesh.voxel_links[k].center_to_center_distance;
        int i = pV1->mesh_index;  // "testing" voxel
        int j = pV2->mesh_index;  // "Target" voxel
        double J_ij;
        
        // k is the Link index!

        double density1 = voxel_population_vectors[i].live_cell_counts[1]/voxel_population_vectors[i].phenotypes_vector[1].max_cells;
        double density2 = voxel_population_vectors[j].live_cell_counts[1]/voxel_population_vectors[j].phenotypes_vector[1].max_cells;


        // FROM 11/2
        
//        int h_ij = tissue_mesh.heaviside_fn(substrate_vectors[j].substrate_quantity[1] - substrate_vectors[i].substrate_quantity[1]);
//
//        int h_ji = tissue_mesh.heaviside_fn(substrate_vectors[i].substrate_quantity[1] - substrate_vectors[j].substrate_quantity[1]);
        //        Variant 1 - "diffusive rho_v" 11-2-17
        
        //        J_ij = -fmax(0,(density1 - density2))*voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])*h_ij
        //        + fmax(0,density2-density1)*voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1])*h_ji;
        
        //        Variant 2 - non-difference rho_v, perhaps more advective 11-29-17
        
//        J_ij = -fmax(0,(density1))*voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])*h_ij
//        + fmax(0,density2)*voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1])*h_ji;
        
        //        std::cout<<voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])<<std::endl;
        // END 11/2
        
        // FROM 11/7 - Why would this not work on a sphere?  - well if the diffusion equations aren't right ... like if you accidently get a reversed gradient ... duh.  Fix that first.

        J_ij = -fmax(0,(substrate_vectors[j].substrate_quantity[1] - substrate_vectors[i].substrate_quantity[1])/(1-0))*voxel_population_vectors[i].update_vascular_creep_rate( substrate_vectors[i].substrate_quantity[1])
        *density1 + fmax(0,(substrate_vectors[i].substrate_quantity[1]-substrate_vectors[j].substrate_quantity[1])/(1-0))
        *voxel_population_vectors[j].update_vascular_creep_rate( substrate_vectors[j].substrate_quantity[1])*density2;
        
        // END 11/7

        // From 11/30/17  - will test later
        
        // Variant 1 - "diffusive rho_v" (like from 11/7) but with angioF scaled by the value in the voxel versus the range of the substart AND taking a density difference instead of being just the denisty
        
        // Varient 2 - "advective rho_v" like 11/7 but with angioF scaled by the value in the voxel versus the range of the substart
        
        //        J_apoptotic_cells[i] = (density2 - density1)/dx;
        //        J_necrotic_cells[i] = (density2 - density1)/dx;

        change_in_live_cell_population[k] = -dt*J_ij*area;

    }

    for( int k=0; k<tissue_mesh.voxel_links.size(); k++)
    {
        Voxel* pV1 = tissue_mesh.voxel_links[k].pVoxel1;  // Is there a better way to store/access there rather than just access then twice in a row????
        Voxel* pV2 = tissue_mesh.voxel_links[k].pVoxel2;
        int i = pV1->mesh_index;
        int j = pV2->mesh_index;

        voxel_population_vectors[i].live_cell_counts[1] = voxel_population_vectors[i].live_cell_counts[1] - change_in_live_cell_population[k];
        voxel_population_vectors[j].live_cell_counts[1] = voxel_population_vectors[j].live_cell_counts[1] + change_in_live_cell_population[k];
    }


    return;
}
    
void Tissue::flux_cell_populations ( double dt )
{
    std::vector<double> change_in_live_cell_population;
    std::vector<double> change_in_apoptotic_cell_population;
    std::vector<double> change_in_necrotic_cell_population;
    
    change_in_live_cell_population.assign(tissue_mesh.voxel_links.size(), 0.0);
    change_in_apoptotic_cell_population.assign(tissue_mesh.voxel_links.size(), 0.0);
    change_in_necrotic_cell_population.assign(tissue_mesh.voxel_links.size(), 0.0);
    
    // Calculate fluxes

    for( int i=0; i<tissue_mesh.voxel_links.size(); i++)
    {
        Voxel* pV1 = tissue_mesh.voxel_links[i].pVoxel1;
        Voxel* pV2 = tissue_mesh.voxel_links[i].pVoxel2;
        double area = tissue_mesh.voxel_links[i].surface_area;
        double dx = tissue_mesh.voxel_links[i].center_to_center_distance;
        int voxel1_index = pV1->mesh_index;
        int voxel2_index = pV2->mesh_index;

        double total_cells1 = (voxel_population_vectors[voxel1_index].live_cell_counts[0]+voxel_population_vectors[voxel1_index].apoptotic_cell_counts[0]+voxel_population_vectors[voxel1_index].necrotic_cell_counts[0]);
        double total_cells2 = (voxel_population_vectors[voxel2_index].live_cell_counts[0]+voxel_population_vectors[voxel2_index].apoptotic_cell_counts[0]+voxel_population_vectors[voxel2_index].necrotic_cell_counts[0]);
        double density1 = total_cells1/voxel_population_vectors[voxel1_index].phenotypes_vector[0].max_cells;
        double density2 = total_cells2/voxel_population_vectors[voxel2_index].phenotypes_vector[0].max_cells;
        
        int h_ij = tissue_mesh.heaviside_fn(total_cells1 - voxel_population_vectors[voxel1_index].phenotypes_vector[0].max_cells*voxel_population_vectors[voxel1_index].phenotypes_vector[0].spatial_mechanical_factor)*
        tissue_mesh.heaviside_fn(density1-density2);

        int h_ji = tissue_mesh.heaviside_fn(total_cells2 - voxel_population_vectors[voxel2_index].phenotypes_vector[0].max_cells*voxel_population_vectors[voxel2_index].phenotypes_vector[0].spatial_mechanical_factor)*
        tissue_mesh.heaviside_fn(density2-density1);
        
        
        change_in_live_cell_population[i] = -dt*((density2 - density1)/dx*(voxel_population_vectors[voxel1_index].phenotypes_vector[0].motility*(h_ij*
                                             voxel_population_vectors[voxel1_index].live_cell_counts[0]/(total_cells1+1E-16) + h_ji*
                                             voxel_population_vectors[voxel2_index].live_cell_counts[0]/(total_cells2+1E-16)))+0)*area;
        
    change_in_apoptotic_cell_population[i] = -dt*((density2 - density1)/dx*(voxel_population_vectors[voxel1_index].phenotypes_vector[0].motility*(h_ij*
                                             voxel_population_vectors[voxel1_index].apoptotic_cell_counts[0]/(total_cells1+1E-16) + h_ji*
                                             voxel_population_vectors[voxel2_index].apoptotic_cell_counts[0]/(total_cells2+1E-16)))+0)*area;

    change_in_necrotic_cell_population[i] = -dt*((density2 - density1)/dx*(voxel_population_vectors[voxel1_index].phenotypes_vector[0].motility*(h_ij*
                                             voxel_population_vectors[voxel1_index].necrotic_cell_counts[0]/(total_cells1+1E-16) + h_ji*
                                             voxel_population_vectors[voxel2_index].necrotic_cell_counts[0]/(total_cells2+1E-16)))+0)*area;
    }

    for( int i=0; i<tissue_mesh.voxel_links.size(); i++)
    {
        Voxel* pV1 = tissue_mesh.voxel_links[i].pVoxel1;  // Is there a better way to store/access there rather than just access then twice in a row????
        Voxel* pV2 = tissue_mesh.voxel_links[i].pVoxel2;
        int voxel1_index = pV1->mesh_index;
        int voxel2_index = pV2->mesh_index;
        
        voxel_population_vectors[voxel1_index].live_cell_counts[0] = voxel_population_vectors[voxel1_index].live_cell_counts[0] - change_in_live_cell_population[i];
        voxel_population_vectors[voxel2_index].live_cell_counts[0] = voxel_population_vectors[voxel2_index].live_cell_counts[0] + change_in_live_cell_population[i];
        voxel_population_vectors[voxel1_index].apoptotic_cell_counts[0] = voxel_population_vectors[voxel1_index].apoptotic_cell_counts[0] - change_in_apoptotic_cell_population[i];
        voxel_population_vectors[voxel2_index].apoptotic_cell_counts[0] = voxel_population_vectors[voxel2_index].apoptotic_cell_counts[0] + change_in_apoptotic_cell_population[i];
        voxel_population_vectors[voxel1_index].necrotic_cell_counts[0] = voxel_population_vectors[voxel1_index].necrotic_cell_counts[0] - change_in_necrotic_cell_population[i];
        voxel_population_vectors[voxel2_index].necrotic_cell_counts[0] = voxel_population_vectors[voxel2_index].necrotic_cell_counts[0] + change_in_necrotic_cell_population[i];
    }
    

    return;
}
    
void Tissue::write_voxel_populations_to_matlab( std::string filename )
{
    int number_of_data_entries = voxel_population_vectors.size();
    int size_of_each_datum = 4; // Live, necrotic, and apopotic cell counts (per voxel) for up to five populations.

    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "voxels_populations" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    // storing data as cols
    for( int i=0; i < number_of_data_entries ; i++ )
    {

        fwrite( (char*) &( voxel_population_vectors[i].live_cell_counts.at(0) ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( voxel_population_vectors[i].apoptotic_cell_counts.at(0) ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( voxel_population_vectors[i].necrotic_cell_counts[0] ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( voxel_population_vectors[i].live_cell_counts.at(1) ) , sizeof(double) , 1 , fp );
      
    }

    fclose( fp );

    return;
}
    
void Tissue::write_substrates_to_matlab( std::string filename )
{
    int number_of_data_entries = substrate_vectors.size();
    int size_of_each_datum = 2;
    
    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "substrates" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    // storing data as cols
    for( int i=0; i < number_of_data_entries ; i++ )
    {
        
        fwrite( (char*) &( substrate_vectors[i].substrate_quantity.at(0) ) , sizeof(double) , 1 , fp );
        fwrite( (char*) &( substrate_vectors[i].substrate_quantity.at(1) ) , sizeof(double) , 1 , fp );

    }

    fclose( fp );

    return;
}
    
void Tissue::write_substrate_properties_to_matlab( std::string filename )
{
    int number_of_data_entries = substrate_properties_vectors.size();
    int size_of_each_datum = 1;
    
    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "substrates" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )
    {
        
        fwrite( (char*) &( substrate_properties_vectors[i].diffusion_coefficients.at(0) ) , sizeof(double) , 1 , fp );

    }

    fclose( fp );

    return;
}

void Tissue::write_all_to_matlab( std::string filename )
{
    int number_of_data_entries = voxel_population_vectors.size();
    int size_of_each_datum = 8;
    
    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "voxels_populations" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )
    {
        
        fwrite( (char*) &( voxel_population_vectors[i].live_cell_counts.at(0) ) , sizeof(double) , 1 , fp );// Live tumor cells
        fwrite( (char*) &( voxel_population_vectors[i].apoptotic_cell_counts.at(0) ) , sizeof(double) , 1 , fp ); //  apoptotic tumor cells
        fwrite( (char*) &( voxel_population_vectors[i].necrotic_cell_counts[0] ) , sizeof(double) , 1 , fp );  //  necrotic tumor cells
        fwrite( (char*) &( voxel_population_vectors[i].live_cell_counts.at(1) ) , sizeof(double) , 1 , fp );  //  Vascular cells
        fwrite( (char*) &( voxel_population_vectors[i].apoptotic_cell_counts.at(1) ) , sizeof(double) , 1 , fp );  //  Vasculars - NOT BEING USED
        fwrite( (char*) &( voxel_population_vectors[i].necrotic_cell_counts[1] ) , sizeof(double) , 1 , fp );   //  Vasular field - NOT BEING USED
        fwrite( (char*) &( substrate_vectors[i].substrate_quantity.at(0) ) , sizeof(double) , 1 , fp );  //  Oxygen
        fwrite( (char*) &( substrate_vectors[i].substrate_quantity.at(1) ) , sizeof(double) , 1 , fp );  //  Anigogenic Factor

    }

    fclose( fp );

    return;
}
    
//struct Tissue_SVG_options_struct Tissue_SVG_options;

void Tissue::write_population_svg (int population_number, double time, std::string filename, double height, double width)
{
    std::ofstream os( filename , std::ios::out );
    
    Write_SVG_start( os, height , width );

    for(int i=0; i<tissue_mesh.voxels.size(); i++)
    {
        
        //double r;
        
        std::vector<std::string> voxel_color;
        
        //r = pow(tissue.tissue_mesh.voxels[i].volume,1/3);
        
        // What will be writing out the links to an image.
        
        //for(int i=0; i<tissue.tissue_mesh.voxel_links.size(); i++)
        //{
        
        //    Write_SVG_line( os , double start_x, double start_y, double end_x , double end_y, 0.25, "rgb(0,0,0)");  // In theory will want to store this information I suppose instead of writing it everytime.  We won't be sprouting voxels and links.
        //}
        
        //                         voxel_color = tissue.blue_to_yellow_coloring( &tissue.voxel_population_vectors[i]);
        voxel_color = population_blue_to_yellow_coloring( population_number, &voxel_population_vectors[i]);
        Write_SVG_circle( os, tissue_mesh.voxels[i].center[0]+750, tissue_mesh.voxels[i].center[1]+750, 19, 0.5, "rgb(0,0,0)" , voxel_color[0]);
        formatted_minutes_to_DDHHMM( time );
//        Write_SVG_text( os, const char* str , double position_x, double position_y, double font_size , const char* color , const char* font);  - This probably needs the structure ...
        //
        //                         std::cout<<tissue.tissue_mesh.voxels[i].center[0]+750<<std::endl;
        //                         std::cout<<tissue.tissue_mesh.voxels[i].center[1]+750<<std::endl;
        
        // Writing out the links.  I suppose ideally there would only happen once and then it would keep being used.  No need to keep figuring out link that don't change.  Also, as written, it will find the links twice.  Probably just test over the voxel links/list of voxel links ... but this is a start ... maybe can use the normal vector or the distance or something ... but for the moment ...
        
        //                         for(int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
        //                         {
        //                             double x_coordinate_i = tissue.tissue_mesh.voxels[i].center[0];
        //                             double y_coordinate_i =tissue.tissue_mesh.voxels[i].center[1];
        //
        //
        //                             for(int j=i; j<tissue.tissue_mesh.voxels.size(); j++)
        //                             {
        //
        //                                 double x_coordinate_j = tissue.tissue_mesh.voxels[j].center[0];
        //                                 double y_coordinate_j = tissue.tissue_mesh.voxels[j].center[1];
        //
        //                                 double distance_sqr = pow(x_coordinate_i - x_coordinate_j,2) + pow(y_coordinate_i-y_coordinate_j,2);
        //                                 double test_distance_sqr = pow(test_distance,2);
        //
        //
        //                                 if((distance_sqr<test_distance_sqr+1E-2)&&(i!=j))
        //                                 {
        //
        //                                     // How will I make the link printing shorter - not from center to center?
        //
        //                                     Write_SVG_line( os , tissue.tissue_mesh.voxels[i].center[0]+750, tissue.tissue_mesh.voxels[i].center[1]+750, tissue.tissue_mesh.voxels[j].center[0] +750, tissue.tissue_mesh.voxels[j].center[1]+750, 2.5, "rgb(0,0,0)");
        //
        //                                     Write_SVG_line( os , tissue.tissue_mesh.voxels[i].center[0]+750, tissue.tissue_mesh.voxels[i].center[1]+750, tissue.tissue_mesh.voxels[j].center[0] +750, tissue.tissue_mesh.voxels[j].center[1]+750, 1.5, "rgb(255,255,255)");
        //                                 }
        //
        //                             }
        //
        //                         }
        
    }
    //std::cout<<voxel_color[0]<<std::endl;
    Write_SVG_end( os );
    
    os.close();
    
    
    return;
}

std::vector<std::string> Tissue::population_blue_to_yellow_coloring( int population_number, BioFVM::Voxel_Population_Vector* pVPV)
{
    
    static std::vector< std::string > output(1, "rgb(0,0,0)");
    
    // Determine live cell density in tissue voxel (live cell count/max cell count)
    
    double cell_density;
    cell_density = pVPV->live_cell_counts[population_number] / pVPV->phenotypes_vector[population_number].max_cells;
    
    // Determine RGB value
    
    double temp1 = 255*cell_density ;
    double temp2 = 255*cell_density;
    double temp3 = (255 * (1-cell_density));
    
    if(temp1< 0)
    {
        std::cout<<"Error in color function - rgb too low"<<std::endl;
        
    }
    
    if(temp1>255)
    {
        std::cout<<"Error in color function - rgb too high"<<std::endl;
        temp1 = 255;
        temp2 = 255;
        temp3 = 0;
    }
    
    if(temp2<0)
    {
        temp1 = 0;
        temp2 = 0;
        temp3 = 255;
    }
    // Print into string output
    
    static char szTempString [128];
    sprintf( szTempString, "rgb(%u,%u,%u)", (int) round( temp1 ), (int) round( temp2 ), (int) round( temp3 ));
    output[0].assign( szTempString );
    //std::cout<<output[0]<<std::endl;
    return output;
}

std::string Tissue::formatted_minutes_to_DDHHMM( double minutes )
{
    static std::string output;
    output.resize( 1024 );
    
    int nMinutes = rint(minutes); // round( minutes );
    // int nDays = (int) floor( (minutes+1e-6) / 1440.0 ); // minutes / 1440
    int nDays = nMinutes / 1440;
    nMinutes -= nDays*1440;
    
    // int nHours = (int) floor( (nMinutes+1e-6) / 60.0 ); // nMinutes / 60;
    int nHours = nMinutes / 60;
    double dMinutes = minutes - 60*( nDays*24 + nHours );
    if( dMinutes < 0 )
    { dMinutes = 0.0; }
    sprintf( (char*) output.c_str(),"%d days, %d hours, and %2.2f minutes", nDays,nHours,dMinutes);
    
    return output ;
}

void Tissue::write_subsrate_svg (int substrate_number, double time, std::string filename, double height, double width)
{
    
    std::ofstream os( filename , std::ios::out );
    
    Write_SVG_start( os, height , width );

    for(int i=0; i<tissue_mesh.voxels.size(); i++)
    {
        
        //double r;
        
        std::vector<std::string> voxel_color;
        
        //r = pow(tissue.tissue_mesh.voxels[i].volume,1/3);
        
        // What will be writing out the links to an image.
        
        //for(int i=0; i<tissue.tissue_mesh.voxel_links.size(); i++)
        //{
        
        //    Write_SVG_line( os , double start_x, double start_y, double end_x , double end_y, 0.25, "rgb(0,0,0)");  // In theory will want to store this information I suppose instead of writing it everytime.  We won't be sprouting voxels and links.
        //}
        
        //                         voxel_color = tissue.blue_to_yellow_coloring( &tissue.voxel_population_vectors[i]);
        voxel_color = substrate_blue_to_yellow_coloring( substrate_number, &substrate_vectors[i]);
//        std::cout<<"voxel colors"<<voxel_color[0]<<std::endl;
        Write_SVG_circle( os, tissue_mesh.voxels[i].center[0]+750, tissue_mesh.voxels[i].center[1]+750, 19, 0.5, "rgb(0,0,0)" , voxel_color[0]);
        formatted_minutes_to_DDHHMM( time );
        //        Write_SVG_text( os, const char* str , double position_x, double position_y, double font_size , const char* color , const char* font);  - This probably needs the structure ...
        //
        //                         std::cout<<tissue.tissue_mesh.voxels[i].center[0]+750<<std::endl;
        //                         std::cout<<tissue.tissue_mesh.voxels[i].center[1]+750<<std::endl;
        
        // Writing out the links.  I suppose ideally there would only happen once and then it would keep being used.  No need to keep figuring out link that don't change.  Also, as written, it will find the links twice.  Probably just test over the voxel links/list of voxel links ... but this is a start ... maybe can use the normal vector or the distance or something ... but for the moment ...
        
        //                         for(int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
        //                         {
        //                             double x_coordinate_i = tissue.tissue_mesh.voxels[i].center[0];
        //                             double y_coordinate_i =tissue.tissue_mesh.voxels[i].center[1];
        //
        //
        //                             for(int j=i; j<tissue.tissue_mesh.voxels.size(); j++)
        //                             {
        //
        //                                 double x_coordinate_j = tissue.tissue_mesh.voxels[j].center[0];
        //                                 double y_coordinate_j = tissue.tissue_mesh.voxels[j].center[1];
        //
        //                                 double distance_sqr = pow(x_coordinate_i - x_coordinate_j,2) + pow(y_coordinate_i-y_coordinate_j,2);
        //                                 double test_distance_sqr = pow(test_distance,2);
        //
        //
        //                                 if((distance_sqr<test_distance_sqr+1E-2)&&(i!=j))
        //                                 {
        //
        //                                     // How will I make the link printing shorter - not from center to center?
        //
        //                                     Write_SVG_line( os , tissue.tissue_mesh.voxels[i].center[0]+750, tissue.tissue_mesh.voxels[i].center[1]+750, tissue.tissue_mesh.voxels[j].center[0] +750, tissue.tissue_mesh.voxels[j].center[1]+750, 2.5, "rgb(0,0,0)");
        //
        //                                     Write_SVG_line( os , tissue.tissue_mesh.voxels[i].center[0]+750, tissue.tissue_mesh.voxels[i].center[1]+750, tissue.tissue_mesh.voxels[j].center[0] +750, tissue.tissue_mesh.voxels[j].center[1]+750, 1.5, "rgb(255,255,255)");
        //                                 }
        //
        //                             }
        //
        //                         }
        
    }
    //std::cout<<voxel_color[0]<<std::endl;
    Write_SVG_end( os );
    
    os.close();
    
    
    return;
    
}

std::vector<std::string> Tissue::substrate_blue_to_yellow_coloring( int substrate_number, BioFVM::Substrate_Vector* pSV)
{

    std::vector< std::string > output(1, "rgb(0,0,0)");
    // Determine live cell density in tissue voxel (live cell count/max cell count)
    
    double substrate_density;
    substrate_density = pSV->substrate_quantity[substrate_number] / 1;
    //std::cout<<pVPV->live_cell_counts[0]<<std::endl;
    //std::cout<<pVPV->phenotypes_vector[0].max_cells<<std::endl;
    //std::cout<<cell_density<<std::endl;
    
    // Determine RGB value
    
    double temp1;
    double temp2;
    double temp3;

//    std::cout<<substrate_density<<std::endl;
//    double temp1 = 255*substrate_density ;  // Not getting the dynamic range I was hoping for ...
//    double temp2 = 255*substrate_density;
//    double temp3 = (255 * (1-substrate_density));

    if(substrate_density<0.001)
    {
         temp1 = 0;  // Not getting the dynamic range I was hoping for ...
         temp2 = 0;
         temp3 = 255;
//        std::cout<<"blue"<<std::endl;
//        std::cout<<substrate_density<<std::endl;
    }

    if(substrate_density>=0.001 && substrate_density<0.5)
    {
         temp1 = 125;  // Not getting the dynamic range I was hoping for ...
         temp2 = 125;
         temp3 = 125;
//        std::cout<<"gray"<<std::endl;
//        std::cout<<substrate_density<<std::endl;
    }

    if(substrate_density>=0.5 && substrate_density<1.0)
    {
         temp1 = 255;  // Not getting the dynamic range I was hoping for ...
         temp2 = 255;
         temp3 = 0;
//        std::cout<<"yellow"<<std::endl;
//        std::cout<<substrate_density<<std::endl;
    }
    
    if(temp1< 0)
    {
        std::cout<<"Error in color function - rgb too low"<<std::endl;
        temp1 = 0;
        temp2 = 0;
        temp3 = 255;
    }

    if(temp1>255)
    {
        std::cout<<"Error in color function - rgb too high"<<std::endl;
        temp1 = 255;
        temp2 = 255;
        temp3 = 0;
    }
    
    if(temp2<0)
    {
        temp1 = 0;
        temp2 = 0;
        temp3 = 255;
    }
    // Print into string output

    static char szTempString [128];
    sprintf( szTempString, "rgb(%u,%u,%u)", (int) round( temp1 ), (int) round( temp2 ), (int) round( temp3 ));

    output[0].assign( szTempString );
    
    return output;
}
    
Cartesian_Tissue::Cartesian_Tissue()
 {
     Cartesian_Mesh mesh;  // Add constructor for using the links!
     mesh.use_voxel_links = "true";
     
     
     return;
 }
    
void Cartesian_Tissue::create_voxel_population_vectors ( void )  //What is happening here?  This doens't make much if any sense.
 {
     
     for(int i = 0; i<mesh.voxels.size(); i++)
     {
         Voxel_Population_Vector population_vector;
         population_vector.tissue_voxel_index = i; // SHould i just make a function call for this instead of using the default constructor?
         population_vector.pMicroenvironment = &chemical_microenvironment;  //
     }
     
     return;
 }

//Tissue_SVG_options_struct Tissue_SVG_options;
//
//std::vector<std::string> Tissue_SVG::black_and_white_coloring( BioFVM::Voxel_Population_Vector* pVPV)
// {
//
//     static std::vector< std::string > output(1, "rgb(0,0,0)");
//
//     // Determine live cell density in tissue voxel (live cell count/max cell count)
//
//     double cell_density;
//     cell_density = pVPV->live_cell_counts[0] / pVPV->phenotypes_vector[0].max_cells;
//
//     // Determine RGB value
//
//     double temp = 255 * cell_density;
//
//     // Print into string output
//
//     static char szTempString [128];
//     springf( szTempString, "rgb(%u,%u,%u)", (int) round( temp ), (int) round( temp ), (int) round( temp ));
//     output[0].assign( szTempString );
//
//     return output:
// }
//
//std::string Tissue_SVG::formatted_minutes_to_DDHHMM( double minutes )
// {
//    static std::string output;
//    output.resize( 1024 );
//
//    int nMinutes = rint(minutes); // round( minutes );
//    // int nDays = (int) floor( (minutes+1e-6) / 1440.0 ); // minutes / 1440
//    int nDays = nMinutes / 1440;
//    nMinutes -= nDays*1440;
//
//    // int nHours = (int) floor( (nMinutes+1e-6) / 60.0 ); // nMinutes / 60;
//    int nHours = nMinutes / 60;
//    double dMinutes = minutes - 60*( nDays*24 + nHours );
//    if( dMinutes < 0 )
//    { dMinutes = 0.0; }
//    sprintf( (char*) output.c_str(),"%d days, %d hours, and %2.2f minutes", nDays,nHours,dMinutes);
//
//    return output ;
// }
    
}

