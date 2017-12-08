#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#ifndef __PhysiCell_vascularization__
#define __PhysiCell_vascularization__

namespace PhysiCell{

class Vascular_Options
{
 private:
 public:
	int vascular_mesh_multiplier; // a multiple of the BioFVM mesh size 
	double vascular_extension_speed; 
	double vascular_birth_rate; 
	
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
	
	Vascular_Densities(); // done 
}; 

class Coarse_Vasculature
{
 private:
 public: 
	Cartesian_Mesh mesh; 
	std::vector<Vascular_Densities> vascular_densities; 
	
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
}; 

void coarse_vasculature_setup( void ); 
void update_coarse_vasculature( double dt ); 

void vascular_supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here );
void vascular_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here );
void vascular_uptake_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here ); 



extern Coarse_Vasculature coarse_vasculature; 

void add_VEGF_to_BioFVM( void ); 
void add_VEGF_to_cells( void ); 

void angiogenesis_bulk_source_function( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination ); 
void angiogenesis_bulk_supply_target_densities_function( Microenvironment* pMicroenvironment, int voxel_index, std::vector<double>* write_destination ); 

 

};

#endif
