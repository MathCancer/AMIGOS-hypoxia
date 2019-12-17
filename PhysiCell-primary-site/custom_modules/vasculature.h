//
//  vasculature.hpp
//  
//
//  Created by John Metzcar on 4/10/18.
//

#ifndef vasculature_h
#define vasculature_h

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

namespace PhysiCell{

        
class Vascular_Options
{
private:
public:
    int vascular_mesh_multiplier; // a multiple of the BioFVM mesh size
    double base_vascular_extension_rate;
    double vascular_birth_rate;
    double vascular_death_rate;
    double max_vascular_density;
    double vascular_proliferation_threshold;
    double vascular_proliferation_saturation;
    double vascular_chemotaxis_threshold;
    double vascular_chemotaxis_saturation;
    double effective_vascular_cutoff_threshold;
    
    double angiogenesis_dt;
    
    double degradation_rate_per_cell;
    
    std::vector<double> blood_substrate_densities;
    std::vector<double> tissue_far_field_substrate_densities;
    
    double blood_oxygen_tension;
    double tissue_far_field_oxygen_tension;
    
    Vascular_Options(); // done
    
    void sync_to_BioFVM( void ); // done
};

extern Vascular_Options default_vascular_options;

class Vascular_Densities
{
private:
public:
    double functional;
    double total;
    double vascular_extension_rate;
    double secretion_rate;
    double target_O2;
    double target_ECM;
	double target_VEGF;
	std::vector<double> target_vector;
	
    Vascular_Densities(); // done
};

class Coarse_Vasculature
{
private:
public:
    Cartesian_Mesh mesh;
    std::vector<Vascular_Densities> vascular_densities;
    
    std::vector<double> VEGF;
    std::vector<double> net_vascular_density_fluxes;
    
    int number_of_voxel_links;
    
    
    Microenvironment* pMicroenvironment;
    
    Vascular_Densities& operator()( std::vector<double> position ); // get densities neariest to (x,y,z)
    Vascular_Densities& operator()( int n ); // done // get densities at vascular index n
    Vascular_Densities& operator()( int i, int j, int k ); // done // get densities at vascular index (i,j,k)
    Vascular_Densities& operator()( Cell* pCell ); // done // get densities at or near the cell's position
    
    // later: add a function that figures out the coarse voxel index based on
    // fine voxel index (n) or fine voxel indices (i,j,k)
    
    
    
    Coarse_Vasculature(); // done
    
    std::vector<double> blood_substrate_densities;
    
    // set to size of BioFVM, add VEGF to TME, check for O2
    //
    void sync_to_BioFVM( void ); // done
    
    void compute_coarse_VEGF( void );
};

void coarse_vasculature_setup( void );

// Adding lines below as temporary mesasure to get around code that is failing in the bulk solver
//std::vector<double> one;
//std::vector<double> zero;
//std::vector<double> bulk_vascular_source_sink_solver_temp1;
//std::vector<double> bulk_vascular_source_sink_solver_temp2;
//std::vector<double> bulk_vascular_source_sink_solver_temp3;
//bool vascular_source_sink_solver_setup_done;

// End ad hoc add. Variables are used in the function below - update_coarse_vasculature

void update_coarse_vasculature( double dt );
void update_vascular_population ( Microenvironment* microenvironment, double dt );
void flux_vascular_density ( Microenvironment* pMicroenvironment, double dt );
void update_vascular_extension_rate ( Microenvironment* pMicroenvironment);

void vascular_supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here );
void vascular_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here );
void vascular_uptake_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here );



extern Coarse_Vasculature coarse_vasculature;
extern Vascular_Densities vascular_densities;

void add_VEGF_to_BioFVM( void );
void add_VEGF_to_cells( void );

void angiogenesis_bulk_source_function( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination );
void angiogenesis_bulk_supply_target_densities_function( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination );
    
void write_vasculature_data_matlab( std::string filename );

    
};


#endif /* vasculature_h */
